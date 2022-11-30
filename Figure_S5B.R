# Skip the import steps if already performed using the code in Figure_S1.R. They reference the same Input_list and AA10 file.
# Import the raw data
Input_list <-
  list.files(path = "./Input/",
             pattern = "*.txt", 
             full.names = TRUE) %>% 
  lapply(read_plus)

Input_list_names <-   
  list.files(path = "./Input/",
             pattern = "*.txt") %>%
  str_replace(pattern = ".txt", 
              replacement = "")

names(Input_list) <- Input_list_names

# Import the primer efficiency file
AA10 <- read_csv("AA10.csv")


# The analysis
OR1A1_1A1primer <- Label_qPCR_data(CqFile = Input_list$`220317_Exp_Cq`, 
                                   Style = "Simple", 
                                   SampleNames = c("MSS-3256_V1.1.2", "MSS-3257_V1.1.2", "MSS-3258_V1.1.2", "MSS-3337_V1.1.2", "MSS-3339_V1.1.2", 
                                                   "MSS-3372_V2.2", "MSS-3374_V2.2", "MSS-3389_V2.2", "MSS-3390_V2.2", 
                                                   "OR1A1_gDNA", "UPwater"), 
                                   GeneNames = c("OR1A1_m", "Acsm4", "Acss2", "Slc25a35"),
                                   Mutation = c("OR1A1_1.1.2", "OR1A1_1.1.2", "OR1A1_1.1.2", "OR1A1_1.1.2", "OR1A1_1.1.2", 
                                                "OR1A1_2.2", "OR1A1_2.2", "OR1A1_2.2", "OR1A1_2.2", "OR1A1_gDNA", "Blank"),
                                   Replicates = 3)

OR1A1_1A1primer_QE <- Quality_Evaluation(OR1A1_1A1primer, 
                                         MeltFile = Input_list$`220317_Exp_Melt`,
                                         RepDiffAllow = 0.6)

OR1A1_1A1primer_QE_Plot <- OR1A1_1A1primer_QE$QE_table %>% Plot_QE(MutationColumn = "Mutation")

Primer_1A1_Analysis <- OR1A1_1A1primer_QE$QE_table %>% 
  filter(Sample %notin% c("OR1A1_gDNA", "UPwater"))

Primer_1A1_Analysis <- Primer_1A1_Analysis %>%
  filter(Pos != "F1") %>%
  select(Gene, Sample, Mutation, Replicate, Cq) %>%
  group_by(Gene, Sample, Mutation) %>%
  summarize(MeanCq = mean(Cq, na.rm = TRUE)) %>%
  left_join(AA10, by = "Gene")

CalibratorV1.1.2 <- Primer_1A1_Analysis %>% 
  filter(Mutation %in% c("OR1A1_2.2")) %>% 
  group_by(Gene) %>% 
  summarize(Calibrator1A1Cq = mean(MeanCq, na.rm = TRUE))

Primer_1A1_Analysis <- Primer_1A1_Analysis %>%
  left_join(CalibratorV1.1.2, by = c("Gene"))

Primer_1A1_Analysis <- AA10 %>% 
  pivot_wider(names_from = Gene, values_from = Primer_Efficiency) %>% 
  select(Acsm4, Acss2, Slc25a35) %>% 
  bind_cols(Primer_1A1_Analysis) %>% 
  select(Gene, Sample, Mutation, Primer_Efficiency, Calibrator1A1Cq, MeanCq, Acsm4, Acss2, Slc25a35)

C1V2 <- CalibratorV1.1.2 %>% filter(Gene == "Acsm4") %>% select(-Gene) %>% rename(Acsm4CalCq = Calibrator1A1Cq)
C2V2 <- CalibratorV1.1.2 %>% filter(Gene == "Acss2") %>% select(-Gene) %>% rename(Acss2CalCq = Calibrator1A1Cq)
C3V2 <- CalibratorV1.1.2 %>% filter(Gene == "Slc25a35") %>% select(-Gene) %>% rename(Slc25a35CalCq = Calibrator1A1Cq)

Primer_1A1_Analysis <- bind_cols(Primer_1A1_Analysis, C1V2, C2V2, C3V2)

CQ1V2 <- Primer_1A1_Analysis %>% filter(Gene == "Acsm4") %>% select(Sample, MeanCq) %>% rename(Acsm4MeanCq = MeanCq)
CQ2V2 <- Primer_1A1_Analysis %>% filter(Gene == "Acss2") %>% select(Sample, MeanCq) %>% rename(Acss2MeanCq = MeanCq)
CQ3V2 <- Primer_1A1_Analysis %>% filter(Gene == "Slc25a35") %>% select(Sample, MeanCq) %>% rename(Slc25a35MeanCq = MeanCq)

Primer_1A1_Analysis <- Primer_1A1_Analysis %>%
  left_join(CQ1V2, by = c("Sample")) %>%
  left_join(CQ2V2, by = c("Sample")) %>%
  left_join(CQ3V2, by = c("Sample"))

SampleSpecificDCt1A1 <- Primer_1A1_Analysis %>%
  rename(Genotype = Mutation) %>%
  mutate(DCtSample = (Acsm4MeanCq*log2(Acsm4) + Acss2MeanCq*log2(Acss2) + Slc25a35MeanCq*log2(Slc25a35))/3 - MeanCq*log2(Primer_Efficiency), 
         DCtCal = (Acsm4CalCq*log2(Acsm4) + Acss2CalCq*log2(Acss2) + Slc25a35CalCq*log2(Slc25a35))/3 - Calibrator1A1Cq*log2(Primer_Efficiency), 
         DDCt = DCtSample -DCtCal) %>% 
  ungroup() %>% 
  group_by(Sample, Gene, Genotype) %>% 
  summarize(DCtSampleMean = mean(DCtSample)) %>% 
  ungroup()

#Genotype-specific for error calculations
GenotypeSpecificDCt1A1 <- SampleSpecificDCt1A1 %>%
  group_by(Gene, Genotype) %>% 
  summarize(DCtGenotypeMean = mean(DCtSampleMean),
            DCtGenotypeSize = length(DCtSampleMean),
            DCtSEM = sd(DCtSampleMean)/sqrt(DCtGenotypeSize)) %>%
  ungroup()

#The mean WT aka V1.1.2 values for each Gene, acting as the correction for each Sample DDCt and each Mutation DDCt; see below
WTGenotypeSpecificDCt1A1 <- GenotypeSpecificDCt1A1 %>%
  filter(Genotype == "OR1A1_2.2") %>%
  rename(DCtWTGenotypeSize = DCtGenotypeSize, 
         DCtWTGenotypeMean = DCtGenotypeMean, 
         DCtWTSEM = DCtSEM) %>%
  select(-Genotype)

#DDCt for each sample; again for plotting purposes... I am more consistent for DDCq and DDCt
DDCt_SampleSpecificDCt1A1 <- SampleSpecificDCt1A1 %>%
  left_join(WTGenotypeSpecificDCt1A1, by = "Gene") %>%
  select(Sample, Gene, Genotype, DCtSampleMean, DCtWTGenotypeMean) %>%
  mutate(DDCt = DCtSampleMean - DCtWTGenotypeMean)

#DDCt for each mutation group
DDCt_GenotypeSpecificDCt1A1 <- GenotypeSpecificDCt1A1 %>%
  left_join(WTGenotypeSpecificDCt1A1, by = "Gene") %>%
  mutate(DDCtMean = DCtGenotypeMean - DCtWTGenotypeMean, 
         DDCtSEMpropagated = sqrt(DCtSEM^2 + DCtWTSEM^2),
         DDCtCI = qt(0.975, df = DCtGenotypeSize + DCtWTGenotypeSize - 2) * DDCtSEMpropagated)

DDCt_1A1 <- DDCt_SampleSpecificDCt1A1 %>%
  left_join(DDCt_GenotypeSpecificDCt1A1, by = c("Gene", "Genotype"))

DDCt_1A1 %>%
  filter(Gene == "OR1A1_m") %>%
  mutate(Genotype = case_when(Genotype == "OR1A1_1.1.2" ~ "1A1-Cherry",
                              TRUE ~ "1A1-GCaMP")) %>%
  ggplot(aes(Genotype, DDCt)) +
  geom_point(aes(color = Genotype), size = 6) + 
  scale_color_manual(values=c("1A1-Cherry" = "lightcoral", "1A1-GCaMP" = "green")) +
  geom_errorbar(aes(ymin = DDCtMean - DDCtCI, ymax = DDCtMean + DDCtCI), width = 0.4, position = position_dodge(0.6), size = 1) +
  labs(y = expression("Log2 Fold Change (-" ~ Delta*Delta ~ "Ct)") ) +
  geom_signif(y_position =c(3.5), comparisons = list(c("1A1-Cherry", "1A1-GCaMP")), annotation = c("**"), color = "black", size = 1, textsize = 10) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.ticks.length = unit(0.15,"inch"),
        axis.ticks = element_line(size = 1),
        axis.text = element_text(size = 15, color = "black"),
        axis.text.x = element_text(vjust = -1),
        axis.title.y = element_text(size = 17),
        legend.position = "none",
        plot.caption = element_text(size = 7, face = "italic"), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5),
        axis.line = element_line(size = 1))

compare_means(DDCt ~ Genotype, DDCt_1A1, method = "t.test", group.by = "Gene")

DDCt_1A1 %>%
  filter(Gene == "OR1A1_m") %>%
  mutate(Genotype = case_when(Genotype == "OR1A1_1.1.2" ~ "1A1-Cherry",
                              TRUE ~ "1A1-GCaMP")) %>%
  ggplot(aes(Genotype, DDCt)) +
  geom_point(aes(color = Genotype), size = 10) + 
  scale_color_manual(values=c("1A1-Cherry" = "lightcoral", "1A1-GCaMP" = "green")) +
  scale_y_continuous(limits = c(-0.5,4))+
  geom_errorbar(aes(ymin = DDCtMean - DDCtCI, ymax = DDCtMean + DDCtCI), width = 0.4, position = position_dodge(0.6), size = 1.5) +
  labs(y = expression("Log2 Fold Change (-" ~ Delta*Delta ~ "Ct)") ) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.ticks.length = unit(0.15,"inch"),
        axis.ticks = element_line(size = 2),
        axis.text = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.caption = element_text(size = 7, face = "italic"), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5),
        axis.line = element_line(size = 2))