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

#----Loading the 5AN1 data-----
OR5AN1_6v15_1 <- Label_qPCR_data(CqFile = Input_list$`220308_Exp_Cq_1`, 
                                 Style = "Custom", 
                                 SampleNames = c(rep(c("MSS-3277", "MSS-3278", "MSS-3279", "MSS-3281", "MSS-4393", "MSS-4394", "MSS-4396", "MSS-4397", "MSS-4398", "MSS-4399", "MSS-4400", "MSS-4401"), 27),
                                                 rep(c("WT_gDNA", "WT_gDNA", "WT_gDNA", "Blank", "Blank", "Blank"), 9)),
                                 GeneNames = c(rep(c("Olfr358", "Olfr690", "Olfr1154", "Acsm4", "Acss2", "Slc25a35", "Olfr596_603", "Olfr390", "Olfr510"), each = 36),
                                               rep(c("Olfr358", "Olfr690", "Olfr1154", "Acsm4", "Acss2", "Slc25a35", "Olfr596_603", "Olfr390", "Olfr510"), each = 6)),
                                 Mutation =    c(rep(c("L15_WT", "L15_5AN1", "L15_5AN1", "L15_5AN1", "L6_WT", "L6_WT", "L6_5AN1", "L6_5AN1", "L6_5AN1", "L6_5AN1", "L6_WT", "L6_WT"), 27),
                                                 rep(c("WT_gDNA", "WT_gDNA", "WT_gDNA", "Blank", "Blank", "Blank"), 9)),
                                 Mouse_line = c(rep(c("L15", "L15", "L15", "L15", "L6", "L6", "L6", "L6", "L6", "L6", "L6", "L6"), 27),
                                                rep("Ctrl", 54)),
                                 Plate = rep("One", times = 378),
                                 ReplicateNames = c(rep(c("One", "Two", "Three"), times = 9, each = 12),
                                                    rep(c("One", "Two", "Three"), times = 18)))
OR5AN1_6v15_QE1 <- Quality_Evaluation(OR5AN1_6v15_1, 
                                      MeltFile = Input_list$`220308_Exp_Melt_1`,
                                      RepDiffAllow = 0.6)
OR5AN1_6v15_Plot1 <- OR5AN1_6v15_QE1$QE_table %>% Plot_QE(MutationColumn = "Mutation")

OR5AN1_6v15_2 <- Label_qPCR_data(CqFile = Input_list$`220308_Exp_Cq_2`, 
                                 Style = "Custom", 
                                 SampleNames = c(rep(c("MSS-3277", "MSS-3278", "MSS-3279", "MSS-3281", "MSS-4393", "MSS-4394", "MSS-4396", "MSS-4397", "MSS-4398", "MSS-4399", "MSS-4400", "MSS-4401"), 27),
                                                 rep(c("WT_gDNA", "WT_gDNA", "WT_gDNA", "Blank", "Blank", "Blank"), 9)),
                                 GeneNames = c(rep(c("Olfr358", "Olfr690", "Olfr1154", "Acsm4", "Acss2", "Slc25a35", "Olfr596_603", "Olfr390", "Olfr510"), each = 36),
                                               rep(c("Olfr358", "Olfr690", "Olfr1154", "Acsm4", "Acss2", "Slc25a35", "Olfr596_603", "Olfr390", "Olfr510"), each = 6)),
                                 Mutation =    c(rep(c("L15_WT", "L15_5AN1", "L15_5AN1", "L15_5AN1", "L6_WT", "L6_WT", "L6_5AN1", "L6_5AN1", "L6_5AN1", "L6_5AN1", "L6_WT", "L6_WT"), 27),
                                                 rep(c("WT_gDNA", "WT_gDNA", "WT_gDNA", "Blank", "Blank", "Blank"), 9)),
                                 Mouse_line = c(rep(c("L15", "L15", "L15", "L15", "L6", "L6", "L6", "L6", "L6", "L6", "L6", "L6"), 27),
                                                rep("Ctrl", 54)),
                                 Plate = rep("Two", times = 378),
                                 ReplicateNames = c(rep(c("One", "Two", "Three"), times = 9, each = 12),
                                                    rep(c("One", "Two", "Three"), times = 18)))
OR5AN1_6v15_QE2 <- Quality_Evaluation(OR5AN1_6v15_2, 
                                      MeltFile = Input_list$`220308_Exp_Melt_2`,
                                      RepDiffAllow = 0.6)
OR5AN1_6v15_Plot2 <- OR5AN1_6v15_QE2$QE_table %>% Plot_QE(MutationColumn = "Mutation")

OR5AN1_6v15_3 <- Label_qPCR_data(CqFile = Input_list$`220308_Exp_Cq_3`, 
                                 Style = "Custom", 
                                 SampleNames = c(rep(c("MSS-3277", "MSS-3278", "MSS-3279", "MSS-3281", "MSS-4393", "MSS-4394", "MSS-4396", "MSS-4397", "MSS-4398", "MSS-4399", "MSS-4400", "MSS-4401"), 27),
                                                 rep(c("WT_gDNA", "WT_gDNA", "WT_gDNA", "Blank", "Blank", "Blank"), 9)),
                                 GeneNames = c(rep(c("Olfr358", "Olfr690", "Olfr1154", "Acsm4", "Acss2", "Slc25a35", "Olfr596_603", "Olfr390", "Olfr510"), each = 36),
                                               rep(c("Olfr358", "Olfr690", "Olfr1154", "Acsm4", "Acss2", "Slc25a35", "Olfr596_603", "Olfr390", "Olfr510"), each = 6)),
                                 Mutation =    c(rep(c("L15_WT", "L15_5AN1", "L15_5AN1", "L15_5AN1", "L6_WT", "L6_WT", "L6_5AN1", "L6_5AN1", "L6_5AN1", "L6_5AN1", "L6_WT", "L6_WT"), 27),
                                                 rep(c("WT_gDNA", "WT_gDNA", "WT_gDNA", "Blank", "Blank", "Blank"), 9)),
                                 Mouse_line = c(rep(c("L15", "L15", "L15", "L15", "L6", "L6", "L6", "L6", "L6", "L6", "L6", "L6"), 27),
                                                rep("Ctrl", 54)),
                                 Plate = rep("Three", times = 378),
                                 ReplicateNames = c(rep(c("One", "Two", "Three"), times = 9, each = 12),
                                                    rep(c("One", "Two", "Three"), times = 18)))
OR5AN1_6v15_QE3 <- Quality_Evaluation(OR5AN1_6v15_3, 
                                      MeltFile = Input_list$`220308_Exp_Melt_3`,
                                      RepDiffAllow = 0.6)
OR5AN1_6v15_Plot3 <- OR5AN1_6v15_QE3$QE_table %>% Plot_QE(MutationColumn = "Mutation")

OR5AN1_6v15_Fullx <- OR5AN1_6v15_QE1$QE_table %>%
  bind_rows(OR5AN1_6v15_QE2$QE_table) %>%
  bind_rows(OR5AN1_6v15_QE3$QE_table)

OR5AN1_6v15_Full <- OR5AN1_6v15_QE1$QE_table %>%
  bind_rows(OR5AN1_6v15_QE2$QE_table) %>%
  bind_rows(OR5AN1_6v15_QE3$QE_table)

OR5AN1_6v15_Full <- OR5AN1_6v15_Full %>%
  filter(FalseCq == FALSE) %>%
  filter(NoPeak == FALSE) %>%
  filter(NoCq == FALSE) %>%
  filter(MaxCq == FALSE) %>%
  filter(NotWithinCqRange == FALSE) %>%
  filter(Mouse_line != "Ctrl") %>%
  filter(Peak == 1)
#---------------------------Contains incomplete code for first three plates----
OR1A1_V1_vs_V2.2_1 <- Label_qPCR_data(CqFile = Input_list$`Exp-Cq-220202`, 
                                      Style = "Simple", 
                                      SampleNames = c("MSS-3184_V2.1", "MSS-3254_V1.1.2", "MSS-3255_V1.1.2", "MSS-3256_V1.1.2", "MSS-3257_V1.1.2", "MSS-3258_V1.1.2",
                                                      "MSS-3337_V1.1.2", "MSS-3338_V1.1.2", "MSS-3339_V1.1.2", "MSS-3340_V1.1.2", 
                                                      "MSS-3372_V2.2", "MSS-3373_V2.2", "MSS-3374_V2.2", "MSS-3388_V2.2", "MSS-3389_V2.2", "MSS-3390_V2.2", 
                                                      "WT_gDNA", "UPwater"), 
                                      GeneNames = c("Olfr358", "Olfr609", "Olfr690", "Acsm4", "Acss2", "Slc25a35"),
                                      Mutation = c("WT_2.1", "WT_1.1.2", "WT_1.1.2", "OR1A1_1.1.2", "OR1A1_1.1.2", "OR1A1_1.1.2", 
                                                   "OR1A1_1.1.2", "WT_1.1.2", "OR1A1_1.1.2", "WT_1.1.2", 
                                                   "OR1A1_2.2", "WT_2.2", "OR1A1_2.2", "WT_2.2", "OR1A1_2.2", "OR1A1_2.2", "WT_gDNA", "Blank"),
                                      Plate = rep("One", times = 18),
                                      Replicates = 3)

OR1A1_V1_vs_V2.2_QE_1 <- Quality_Evaluation(OR1A1_V1_vs_V2.2_1, 
                                            MeltFile = Input_list$`Exp-Melt-220202`,
                                            RepDiffAllow = 0.6)

OR1A1_V1_vs_V2.2_QE_Plot_1 <- OR1A1_V1_vs_V2.2_QE_1$QE_table %>% Plot_QE(MutationColumn = "Mutation")

#---------------------------Contains incomplete code for second of three plates----

OR1A1_V1_vs_V2.2_2 <- Label_qPCR_data(CqFile = Input_list$`Exp-Cq-220203_1`, 
                                      Style = "Simple", 
                                      SampleNames = c("MSS-3184_V2.1", "MSS-3254_V1.1.2", "MSS-3255_V1.1.2", "MSS-3256_V1.1.2", "MSS-3257_V1.1.2", "MSS-3258_V1.1.2",
                                                      "MSS-3337_V1.1.2", "MSS-3338_V1.1.2", "MSS-3339_V1.1.2", "MSS-3340_V1.1.2", 
                                                      "MSS-3372_V2.2", "MSS-3373_V2.2", "MSS-3374_V2.2", "MSS-3388_V2.2", "MSS-3389_V2.2", "MSS-3390_V2.2", 
                                                      "WT_gDNA", "UPwater"), 
                                      GeneNames = c("Olfr358", "Olfr609", "Olfr690", "Acsm4", "Acss2", "Slc25a35"),
                                      Mutation = c("WT_2.1", "WT_1.1.2", "WT_1.1.2", "OR1A1_1.1.2", "OR1A1_1.1.2", "OR1A1_1.1.2", 
                                                   "OR1A1_1.1.2", "WT_1.1.2", "OR1A1_1.1.2", "WT_1.1.2", 
                                                   "OR1A1_2.2", "WT_2.2", "OR1A1_2.2", "WT_2.2", "OR1A1_2.2", "OR1A1_2.2", "WT_gDNA", "Blank"),
                                      Plate = rep("Two", times = 18),
                                      Replicates = 3)

OR1A1_V1_vs_V2.2_QE_2 <- Quality_Evaluation(OR1A1_V1_vs_V2.2_2, 
                                            MeltFile = Input_list$`Exp-Melt-220203_1`,
                                            RepDiffAllow = 0.6)

OR1A1_V1_vs_V2.2_QE_Plot_2 <- OR1A1_V1_vs_V2.2_QE_2$QE_table %>% Plot_QE(MutationColumn = "Mutation")

#---------------------------Contains incomplete code for third of three plates----

OR1A1_V1_vs_V2.2_3 <- Label_qPCR_data(CqFile = Input_list$`Exp-Cq-220203_2`, 
                                      Style = "Simple", 
                                      SampleNames = c("MSS-3184_V2.1", "MSS-3254_V1.1.2", "MSS-3255_V1.1.2", "MSS-3256_V1.1.2", "MSS-3257_V1.1.2", "MSS-3258_V1.1.2",
                                                      "MSS-3337_V1.1.2", "MSS-3338_V1.1.2", "MSS-3339_V1.1.2", "MSS-3340_V1.1.2", 
                                                      "MSS-3372_V2.2", "MSS-3373_V2.2", "MSS-3374_V2.2", "MSS-3388_V2.2", "MSS-3389_V2.2", "MSS-3390_V2.2", 
                                                      "WT_gDNA", "UPwater"), 
                                      GeneNames = c("Olfr358", "Olfr609", "Olfr690", "Acsm4", "Acss2", "Slc25a35"),
                                      Mutation = c("WT_2.1", "WT_1.1.2", "WT_1.1.2", "OR1A1_1.1.2", "OR1A1_1.1.2", "OR1A1_1.1.2", 
                                                   "OR1A1_1.1.2", "WT_1.1.2", "OR1A1_1.1.2", "WT_1.1.2", 
                                                   "OR1A1_2.2", "WT_2.2", "OR1A1_2.2", "WT_2.2", "OR1A1_2.2", "OR1A1_2.2", "WT_gDNA", "Blank"),
                                      Plate = rep("Three", times = 18),
                                      Replicates = 3)

OR1A1_V1_vs_V2.2_QE_3 <- Quality_Evaluation(OR1A1_V1_vs_V2.2_3, 
                                            MeltFile = Input_list$`Exp-Melt-220203_2`,
                                            RepDiffAllow = 0.6)

OR1A1_V1_vs_V2.2_QE_Plot_3 <- OR1A1_V1_vs_V2.2_QE_3$QE_table %>% Plot_QE(MutationColumn = "Mutation")

#----This is for the second set of three plates----
OR1A1_V1_vs_V2.2_4 <- Label_qPCR_data(CqFile = Input_list$`Exp-Cq-220218_1`, 
                                      Style = "Simple", 
                                      SampleNames = c("MSS-3184_V2.1", "MSS-3254_V1.1.2", "MSS-3255_V1.1.2", "MSS-3256_V1.1.2", "MSS-3257_V1.1.2", "MSS-3258_V1.1.2",
                                                      "MSS-3337_V1.1.2", "MSS-3338_V1.1.2", "MSS-3339_V1.1.2", "MSS-3340_V1.1.2", 
                                                      "MSS-3372_V2.2", "MSS-3373_V2.2", "MSS-3374_V2.2", "MSS-3388_V2.2", "MSS-3389_V2.2", "MSS-3390_V2.2", 
                                                      "WT_gDNA", "UPwater"), 
                                      GeneNames = c("Olfr1154", "Acsm4", "Acss2", "Slc25a35", "Olfr596_603", "Olfr390", "Olfr510"),
                                      Mutation = c("WT_2.1", "WT_1.1.2", "WT_1.1.2", "OR1A1_1.1.2", "OR1A1_1.1.2", "OR1A1_1.1.2", 
                                                   "OR1A1_1.1.2", "WT_1.1.2", "OR1A1_1.1.2", "WT_1.1.2", 
                                                   "OR1A1_2.2", "WT_2.2", "OR1A1_2.2", "WT_2.2", "OR1A1_2.2", "OR1A1_2.2", "WT_gDNA", "Blank"),
                                      Plate = rep("Four", times = 18),
                                      Replicates = 3)

OR1A1_V1_vs_V2.2_QE_4 <- Quality_Evaluation(OR1A1_V1_vs_V2.2_4, 
                                            MeltFile = Input_list$`Exp-Melt-220218_1`,
                                            RepDiffAllow = 0.6)

OR1A1_V1_vs_V2.2_QE_Plot_4 <- OR1A1_V1_vs_V2.2_QE_4$QE_table %>% Plot_QE(MutationColumn = "Mutation")

#----Plate 2 of 3 for the second set----
OR1A1_V1_vs_V2.2_5 <- Label_qPCR_data(CqFile = Input_list$`Exp-Cq-220218_2`, 
                                      Style = "Simple", 
                                      SampleNames = c("MSS-3184_V2.1", "MSS-3254_V1.1.2", "MSS-3255_V1.1.2", "MSS-3256_V1.1.2", "MSS-3257_V1.1.2", "MSS-3258_V1.1.2",
                                                      "MSS-3337_V1.1.2", "MSS-3338_V1.1.2", "MSS-3339_V1.1.2", "MSS-3340_V1.1.2", 
                                                      "MSS-3372_V2.2", "MSS-3373_V2.2", "MSS-3374_V2.2", "MSS-3388_V2.2", "MSS-3389_V2.2", "MSS-3390_V2.2", 
                                                      "WT_gDNA", "UPwater"), 
                                      GeneNames = c("Olfr1154", "Acsm4", "Acss2", "Slc25a35", "Olfr596_603", "Olfr390", "Olfr510"),
                                      Mutation = c("WT_2.1", "WT_1.1.2", "WT_1.1.2", "OR1A1_1.1.2", "OR1A1_1.1.2", "OR1A1_1.1.2", 
                                                   "OR1A1_1.1.2", "WT_1.1.2", "OR1A1_1.1.2", "WT_1.1.2", 
                                                   "OR1A1_2.2", "WT_2.2", "OR1A1_2.2", "WT_2.2", "OR1A1_2.2", "OR1A1_2.2", "WT_gDNA", "Blank"),
                                      Plate = rep("Five", times = 18),
                                      Replicates = 3)

OR1A1_V1_vs_V2.2_QE_5 <- Quality_Evaluation(OR1A1_V1_vs_V2.2_5, 
                                            MeltFile = Input_list$`Exp-Melt-220218_2`,
                                            RepDiffAllow = 0.6)

OR1A1_V1_vs_V2.2_QE_Plot_5 <- OR1A1_V1_vs_V2.2_QE_5$QE_table %>% Plot_QE(MutationColumn = "Mutation")

#-----Plate 3 of 3----
OR1A1_V1_vs_V2.2_6 <- Label_qPCR_data(CqFile = Input_list$`Exp-Cq-220218_3`, 
                                      Style = "Simple", 
                                      SampleNames = c("MSS-3184_V2.1", "MSS-3254_V1.1.2", "MSS-3255_V1.1.2", "MSS-3256_V1.1.2", "MSS-3257_V1.1.2", "MSS-3258_V1.1.2",
                                                      "MSS-3337_V1.1.2", "MSS-3338_V1.1.2", "MSS-3339_V1.1.2", "MSS-3340_V1.1.2", 
                                                      "MSS-3372_V2.2", "MSS-3373_V2.2", "MSS-3374_V2.2", "MSS-3388_V2.2", "MSS-3389_V2.2", "MSS-3390_V2.2", 
                                                      "WT_gDNA", "UPwater"), 
                                      GeneNames = c("Olfr1154", "Acsm4", "Acss2", "Slc25a35", "Olfr596_603", "Olfr390", "Olfr510"),
                                      Mutation = c("WT_2.1", "WT_1.1.2", "WT_1.1.2", "OR1A1_1.1.2", "OR1A1_1.1.2", "OR1A1_1.1.2", 
                                                   "OR1A1_1.1.2", "WT_1.1.2", "OR1A1_1.1.2", "WT_1.1.2", 
                                                   "OR1A1_2.2", "WT_2.2", "OR1A1_2.2", "WT_2.2", "OR1A1_2.2", "OR1A1_2.2", "WT_gDNA", "Blank"),
                                      Plate = rep("Six", times = 18),
                                      Replicates = 3)

OR1A1_V1_vs_V2.2_QE_6 <- Quality_Evaluation(OR1A1_V1_vs_V2.2_6, 
                                            MeltFile = Input_list$`Exp-Melt-220218_3`,
                                            RepDiffAllow = 0.6)

OR1A1_V1_vs_V2.2_QE_Plot_6 <- OR1A1_V1_vs_V2.2_QE_6$QE_table %>% Plot_QE(MutationColumn = "Mutation")

#----Combining all 6 QE full plates prior to filtering the set down----
Total_OR1A1_V1_vs_V2.2_QE <- OR1A1_V1_vs_V2.2_QE_1$QE_table %>%
  bind_rows(OR1A1_V1_vs_V2.2_QE_2$QE_table) %>%
  bind_rows(OR1A1_V1_vs_V2.2_QE_3$QE_table) %>%
  bind_rows(OR1A1_V1_vs_V2.2_QE_4$QE_table) %>%
  bind_rows(OR1A1_V1_vs_V2.2_QE_5$QE_table) %>%
  bind_rows(OR1A1_V1_vs_V2.2_QE_6$QE_table)

#Stepwise filtering of the results
Final_OR1A1 <- Total_OR1A1_V1_vs_V2.2_QE %>%
  filter(FalseCq == FALSE) %>%
  filter(NoPeak == FALSE) %>%
  filter(NoCq == FALSE) %>%
  filter(MaxCq == FALSE) %>%
  filter(Cq < 29) %>%
  select(-MaxCq, -NoCq, -NoPeak, -FalseCq, -NotWithinCqRange) %>%
  filter(Height > 0.1) %>%
  filter(Sample != "WT_gDNA") %>%
  filter(Sample != "UPwater") %>%
  filter(Peak == 1)

#----Testing a combined plate?----
OR1A1x5AN1 <- Final_OR1A1 %>%
  filter(Mutation %in% c("WT_2.2", "WT_1.1.2", "OR1A1_2.2")) %>%
  filter(Gene %notin% "Olfr609") %>%
  select(Plate, Gene, Cq, Sample, Replicate, Mutation) %>%
  group_by(Plate, Gene, Sample, Mutation) %>%
  summarize(MeanCq = mean(Cq, na.rm = TRUE)) %>%
  left_join(AA10, by = "Gene") %>%
  ungroup()

OR5AN1x1A1 <- OR5AN1_6v15_Full %>%
  filter(Mutation %in% c("L15_WT", "L15_5AN1", "L6_WT")) %>%
  select(Plate, Gene, Cq, Sample, Replicate, Mutation) %>%
  mutate(Plate = case_when(Plate == "One" ~ "A",
                           Plate == "Two" ~ "B",
                           Plate == "Three" ~ "C")) %>%
  group_by(Plate, Gene, Sample, Mutation) %>%
  summarize(MeanCq = mean(Cq, na.rm = TRUE)) %>%
  left_join(AA10, by = "Gene") %>%
  ungroup()

ORx1A1_5AN1 <- OR1A1x5AN1 %>%
  bind_rows(OR5AN1x1A1)

Cal_ORx <- ORx1A1_5AN1 %>% 
  filter(Mutation %in% c("WT_2.2", "WT_1.1.2", "L15_WT", "L6_WT")) %>% 
  group_by(Plate, Gene) %>% 
  summarize(CalibratorCq = mean(MeanCq, na.rm = TRUE))

ORx1A1_5AN1 <- ORx1A1_5AN1 %>%
  left_join(Cal_ORx, by = c("Plate", "Gene"))

ORx1A1_5AN1 <- AA10 %>% 
  pivot_wider(names_from = Gene, values_from = Primer_Efficiency) %>% 
  select(Acsm4, Acss2, Slc25a35) %>% 
  bind_cols(ORx1A1_5AN1) %>% 
  select(Plate, Gene, Sample, Mutation, Primer_Efficiency, CalibratorCq, MeanCq, Acsm4, Acss2, Slc25a35)

OlfrC1x <- Cal_ORx %>% filter(Gene == "Acsm4") %>% ungroup() %>% select(-Gene) %>% rename(Acsm4CalCq = CalibratorCq)
OlfrC2x <- Cal_ORx %>% filter(Gene == "Acss2") %>% ungroup() %>% select(-Gene) %>% rename(Acss2CalCq = CalibratorCq)
OlfrC3x <- Cal_ORx %>% filter(Gene == "Slc25a35") %>% ungroup() %>% select(-Gene) %>% rename(Slc25a35CalCq = CalibratorCq)

ORx1A1_5AN1 <- ORx1A1_5AN1 %>%
  left_join(OlfrC1x, by = c("Plate")) %>%
  left_join(OlfrC2x, by = c("Plate")) %>%
  left_join(OlfrC3x, by = c("Plate"))

OlfrCQ1x <- ORx1A1_5AN1 %>% filter(Gene == "Acsm4") %>% ungroup() %>% select(Plate, Sample, MeanCq) %>% rename(Acsm4MeanCq = MeanCq)
OlfrCQ2x <- ORx1A1_5AN1 %>% filter(Gene == "Acss2") %>% ungroup() %>% select(Plate, Sample, MeanCq) %>% rename(Acss2MeanCq = MeanCq)
OlfrCQ3x <- ORx1A1_5AN1 %>% filter(Gene == "Slc25a35") %>% ungroup() %>% select(Plate, Sample, MeanCq) %>% rename(Slc25a35MeanCq = MeanCq)

ORx1A1_5AN1 <- ORx1A1_5AN1 %>%
  left_join(OlfrCQ1x, by = c("Plate", "Sample")) %>%
  left_join(OlfrCQ2x, by = c("Plate", "Sample")) %>%
  left_join(OlfrCQ3x, by = c("Plate", "Sample"))

#----New DDCq method again----
ORx1A1_5AN1 <- ORx1A1_5AN1 %>%
  mutate(Mutation2 = ifelse(Mutation %in% c("L15_WT", "L6_WT", "WT_1.1.2", "WT_2.2"), "WT", Mutation))

#Sample-specific, for plotting purposes
Final_GCaMP_V2_SampleSpecificDCt <- ORx1A1_5AN1 %>% 
  mutate(DCqSample = (Acsm4MeanCq*log2(Acsm4) + Acss2MeanCq*log2(Acss2) + Slc25a35MeanCq*log2(Slc25a35))/3 - MeanCq*log2(Primer_Efficiency), 
         DCqCal = (Acsm4CalCq*log2(Acsm4) + Acss2CalCq*log2(Acss2) + Slc25a35CalCq*log2(Slc25a35))/3 - CalibratorCq*log2(Primer_Efficiency), 
         DDCq = DCqSample -DCqCal) %>% 
  ungroup() %>% 
  group_by(Sample, Gene, Mutation2) %>% 
  summarize(DCqSampleMean = mean(DCqSample), 
            DCqSampleSEM = sd(DCqSample)/sqrt(length(DCqSample)),
            DCqSEMsq =  DCqSampleSEM^2) %>% 
  ungroup()

#Mutation-specific for error calculations
Final_GCaMP_V2_MutationSpecificDCt <- ORx1A1_5AN1 %>% 
  mutate(DCqSample = (Acsm4MeanCq*log2(Acsm4) + Acss2MeanCq*log2(Acss2) + Slc25a35MeanCq*log2(Slc25a35))/3 - MeanCq*log2(Primer_Efficiency), 
         DCqCal = (Acsm4CalCq*log2(Acsm4) + Acss2CalCq*log2(Acss2) + Slc25a35CalCq*log2(Slc25a35))/3 - CalibratorCq*log2(Primer_Efficiency), 
         DDCq = DCqSample -DCqCal) %>% 
  ungroup() %>% 
  group_by(Sample, Gene, Mutation2) %>% 
  summarize(DCqSampleMean = mean(DCqSample), 
            DCqSampleSEM = sd(DCqSample)/sqrt(length(DCqSample)),
            DCqSEMsq =  DCqSampleSEM^2) %>% 
  ungroup() %>% 
  group_by(Gene, Mutation2) %>% 
  summarize(DCqMutationMean = mean(DCqSampleMean),
            DCqMutationSize = length(DCqSampleMean),
            DCqSEMnoprop = sd(DCqSampleMean)/sqrt(DCqMutationSize), 
            DCqSEMprop = sqrt(sum(DCqSEMsq))) %>%
  ungroup()

#The mean WT values for each Gene, acting as the correction for each Sample DDCt and each Mutation DDCt; see below
Final_GCaMP_V2_WT_DCt <- Final_GCaMP_V2_MutationSpecificDCt %>%
  filter(Mutation2 == "WT") %>%
  rename(DCqWTSize = DCqMutationSize, 
         DCqWTMean = DCqMutationMean, 
         DCqWTSEMnoprop = DCqSEMnoprop, 
         DCqWTSEMprop = DCqSEMprop) %>%
  select(-Mutation2)

#DDCt for each sample; again for plotting purposes... I am not consistent for DDCq and DDCt
Final_GCaMP_V2_Sample_MeanWT_DDCt <- Final_GCaMP_V2_SampleSpecificDCt %>%
  left_join(Final_GCaMP_V2_WT_DCt, by = "Gene") %>%
  select(Sample, Gene, Mutation2, DCqSampleMean, DCqWTMean) %>%
  mutate(DDCq = DCqSampleMean - DCqWTMean)

#DDCt for each mutation group
Final_GCaMP_V2_Mutation_MeanWT_DDCt <- Final_GCaMP_V2_MutationSpecificDCt %>%
  left_join(Final_GCaMP_V2_WT_DCt, by = "Gene") %>%
  mutate(DDCqMean = DCqMutationMean - DCqWTMean, 
         DDCqSEM1prop = sqrt(DCqSEMnoprop^2 + DCqWTSEMnoprop^2),
         DDCqSEM2prop = sqrt(DCqSEMprop^2 + DCqWTSEMprop^2),
         DDCqCI1prop = qt(0.975, df = DCqMutationSize + DCqWTSize - 2) * DDCqSEM1prop,
         DDCqCI2prop = qt(0.975, df = DCqMutationSize + DCqWTSize - 2) * DDCqSEM2prop)

#DDCt for each mutation group and error bars, without grouping for genes, limited to just Olfrs
Final_GCaMP_V2_MutationNoGeneDCt <- ORx1A1_5AN1 %>% 
  mutate(DCqSample = (Acsm4MeanCq*log2(Acsm4) + Acss2MeanCq*log2(Acss2) + Slc25a35MeanCq*log2(Slc25a35))/3 - MeanCq*log2(Primer_Efficiency), 
         DCqCal = (Acsm4CalCq*log2(Acsm4) + Acss2CalCq*log2(Acss2) + Slc25a35CalCq*log2(Slc25a35))/3 - CalibratorCq*log2(Primer_Efficiency), 
         DDCq = DCqSample -DCqCal) %>% 
  ungroup() %>% 
  group_by(Sample, Gene, Mutation2) %>% 
  summarize(DCqSampleMean = mean(DCqSample), 
            DCqSampleSEM = sd(DCqSample)/sqrt(length(DCqSample)),
            DCqSEMsq =  DCqSampleSEM^2) %>% 
  ungroup() %>% 
  filter(Gene %notin% c("Acsm4", "Acss2", "Slc25a35")) %>%
  group_by(Mutation2) %>% 
  summarize(DCqMutationMean = mean(DCqSampleMean),
            DCqMutationSize = length(DCqSampleMean),
            DCqSEMnoprop = sd(DCqSampleMean)/sqrt(DCqMutationSize), 
            DCqSEMprop = sqrt(sum(DCqSEMsq))) %>%
  ungroup()

#The No Gene group DDCt needs its own WT Mean and associated values Error values, calculated here
Final_GCaMP_V2_WTNoGeneDCt <- Final_GCaMP_V2_MutationNoGeneDCt %>%
  filter(Mutation2 == "WT") %>%
  rename(DCqWTSize = DCqMutationSize, 
         DCqWTMean = DCqMutationMean, 
         DCqWTSEMnoprop = DCqSEMnoprop, 
         DCqWTSEMprop = DCqSEMprop) %>%
  select(-Mutation2)

Final_GCaMP_V2_Mutation_WTmean_NoGeneDDCt <- bind_cols(Final_GCaMP_V2_MutationNoGeneDCt, Final_GCaMP_V2_WTNoGeneDCt) %>%
  mutate(DDCqMean = DCqMutationMean - DCqWTMean, 
         DDCqSEM1prop = sqrt(DCqSEMnoprop^2 + DCqWTSEMnoprop^2),
         DDCqSEM2prop = sqrt(DCqSEMprop^2 + DCqWTSEMprop^2),
         DDCqCI1prop = qt(0.975, df = DCqMutationSize + DCqWTSize - 2) * DDCqSEM1prop,
         DDCqCI2prop = qt(0.975, df = DCqMutationSize + DCqWTSize - 2) * DDCqSEM2prop)


#----Plotting for GCaMP----
#The final plotting tables need to combine the Sample table for the points and the Error values, either for each Mutation/Gene group or each Mutation group
#This one for the split Gene Mutation plots, following below
Sample_Mutation_Gene_GCaMP_DDCt <- Final_GCaMP_V2_Sample_MeanWT_DDCt %>%
  left_join(Final_GCaMP_V2_Mutation_MeanWT_DDCt, by = c("Gene", "Mutation2")) %>%
  rename(Genotype = Mutation2) %>%
  mutate(Gene = factor(Gene, levels = c("Acsm4", "Acss2", "Slc25a35", "Olfr358", "Olfr390", "Olfr510", "Olfr596_603", "Olfr609", "Olfr690", "Olfr1154")),
         Genotype = case_when(Genotype == "L15_5AN1" ~ "OR5AN1 Line 15",
                              Genotype == "OR1A1_2.2" ~ "OR1A1 V2.2",
                              TRUE ~ "WT"))

#Full plot with single-prop CI
L2FCGCaMP_1 <- Sample_Mutation_Gene_GCaMP_DDCt %>%
  ggplot(aes(Gene, DDCq, color = Genotype)) +
  geom_point(position = position_dodge(width = 0.6)) + 
  geom_errorbar(aes(ymin = DDCqMean - DDCqCI1prop, ymax = DDCqMean + DDCqCI1prop), width = 0.4, position = position_dodge(0.6)) +
  labs(title = "Log2 Fold Change in Ct values for Reference and Olfr mRNA in OR1A1 V2.2 and OR5AN1 Line 15", 
       subtitle = "Plotted L2FC of each Sample with the 95% confidence interval based on the mean and single-propagated error of each Genotype", 
       caption = "DDCt calculated using corrected primer efficiencies and reference primers Acsm4, Acss2, and Slc25a35. Each Sample and Genotype DDCt calibrated using the appropriate mean WT DCt. Ct and DCt in triplicate.", 
       y = "Log2 Fold Change (DDCt)") +
  theme_bw() +
  theme(axis.title.x = element_blank(), legend.position = c(0.07, 0.15),
        legend.background = element_rect(fill = "white", color = "black"), 
        plot.caption = element_text(size = 7, face = "italic"), 
        plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5))

#First statistical tests
TableGCaMP <- compare_means(DDCq ~ Genotype, Sample_Mutation_Gene_GCaMP_DDCt, method = "t.test", group.by = c("Gene")) %>% 
  filter(Gene %notin% c("Acsm4", "Acss2", "Slc25a35")) %>%
  select(Gene, group1, group2, p) %>%
  arrange(p) %>%
  mutate(row = row_number(),
         HolmAdjustment = 0.05/(18 - row + 1), 
         PassHA = p < HolmAdjustment, 
         HolmCorrection = p * (18 - row + 1),
         PassHC = HolmCorrection < 0.05, 
         Pass_DunnSidakCorrection = p < 1 - (1 - 0.05)^(1/18)) %>% 
  select(Gene, group1, group2, p, HolmCorrection, PassHC, Pass_DunnSidakCorrection) %>% 
  rename(p_value = p, Pass_HolmCorrection = PassHC, Group1 = group1, Group2 = group2) %>%
  arrange(Group1, Group2)

TableGrobGCaMP <- tableGrob(TableGCaMP, rows=NULL)

TextGrobGCaMP <- text_grob(paste("Eighteen independent t.tests were performed. With alpha = 0.05, Holm and Dunn-Šidák corrections were calculated. ", 
                                 "Dunn-Šidák corrected alpha evaluated to 0.00285. The 5AN1 is significantly different from WT for 4/6 Olfr mRNA values.", 
                                 "The 1A1 line is significantly different from WT for 3/6 Olfr mRNA values. The non-WT lines differ statistically for Olfr690.", sep = "\n"), 
                           face = "italic")

grid.arrange(L2FCGCaMP_1, TableGrobGCaMP, TextGrobGCaMP, heights = c(0.7, 0.5, 0.1), as.table = TRUE)

#Full plot with single-prop CI plus combined graph, No References visible, but otherwise the same as above.
L2FCGCaMP_2 <- Sample_Mutation_Gene_GCaMP_DDCt %>%
  filter(Gene %notin% c("Acsm4", "Acss2", "Slc25a35")) %>%
  mutate(Genotype = case_when(Genotype == "OR5AN1 Line 15" ~ "OR5AN1", 
                              Genotype == "OR1A1 V2.2" ~ "OR1A1", 
                              TRUE ~ "WT")) %>%
  mutate(Genotype = factor(Genotype, levels = c("WT", "OR5AN1", "OR1A1"))) %>%
  ggplot(aes(Gene, DDCq, color = Genotype, group = Genotype)) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  scale_color_manual(values=c("WT" = "green", "OR5AN1" = "blue", "OR1A1" = "red")) +
  geom_errorbar(aes(ymin = DDCqMean - DDCqCI1prop, ymax = DDCqMean + DDCqCI1prop), width = 0.4, position = position_dodge(0.6), color = "black") +
  labs(title = "Log2 Fold Change in Ct values for Olfr mRNA in OR5AN1 and OR1A1", 
       subtitle = "Plotted L2FC of each Sample with the 95% confidence interval based on the mean and single-propagated error of each Genotype", 
       caption = "DDCt calculated using corrected primer efficiencies and reference primers Acsm4, Acss2, and Slc25a35. \nEach Sample and Genotype DDCt calibrated using the appropriate mean WT DCt. Ct and DCt in triplicate.", 
       y = "Log2 Fold Change (DDCt)") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(face="bold", size = 10),
        legend.position = c(0.92, 0.15),
        legend.background = element_rect(fill = "white", color = "black"), 
        plot.caption = element_text(size = 7, face = "italic"), 
        plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5))

#Removing some text
IntermediateData1 <- Sample_Mutation_Gene_GCaMP_DDCt %>%
  filter(Gene %notin% c("Acsm4", "Acss2", "Slc25a35")) %>%
  mutate(Genotype = case_when(Genotype == "OR5AN1 Line 15" ~ "OR5AN1", 
                              Genotype == "OR1A1 V2.2" ~ "OR1A1", 
                              TRUE ~ "WT")) %>%
  mutate(Genotype = factor(Genotype, levels = c("WT", "OR5AN1", "OR1A1"))) %>%
  mutate(Class = case_when(Gene %in% c("Olfr596_603", "Olfr690") ~ "Class I", 
                           TRUE ~ "Class II")) %>%
  left_join(DVI2, by = c("Gene" = "Symbol")) %>%
  mutate(DVI = ifelse(is.na(DVI), 1.05, DVI)) %>%
  mutate(Gene = factor(Gene, levels = c("Olfr596_603", "Olfr690", "Olfr510", "Olfr1154", "Olfr358", "Olfr390")))

L2FCGCaMP_3 <-
  ggplot(data = IntermediateData1) +
  geom_point(aes(Gene, DDCq, color = Genotype, group = Genotype), position = position_dodge(width = 0.6), size = 3) +
  scale_color_manual(values=c("WT" = "green3", "OR5AN1" = "blue", "OR1A1" = "red")) +
  geom_errorbar(aes(x = Gene, group = Genotype, ymin = DDCqMean - DDCqCI1prop, ymax = DDCqMean + DDCqCI1prop), 
                width = 0.4, position = position_dodge(0.6), color = "black") +
  labs(title = "Log2 Fold Change in Ct values for Olfr mRNA in OR5AN1 and OR1A1", 
       y = expression("Log2 Fold Change (-" ~ Delta*Delta ~ "Ct)") ) +
  theme_bw() +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text = element_text(face="bold", size = 15),
        axis.title.y = element_text(size = 17),
        legend.position = c(0.92, 0.18),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.text = element_text(face="bold"),
        legend.title = element_blank()) +
  # geom_rect(data = rects, aes(ymin= -3, ymax = 1, xmin=xstart, xmax=xend, fill = ), alpha = 0.25) +
  #geom_rect(data = rects, aes(ymin= -3, ymax = 1, xmin=0.5, xmax=2.5), alpha = 0.05, fill = "blue") +
  #geom_rect(data = rects, aes(ymin= -3, ymax = 1, xmin=2.5, xmax=6.5), alpha = 0.05, fill = "red") +
  geom_vline(xintercept = 2.5, linetype="dotted", color = "grey40", size=1.5) +
  annotate("text", label = "Class I", x = 1.5, y = 0.90, size = 8, colour = "blue", fontface = 2) +
  annotate("text", label = "Class II", x = 4.5, y = 0.90, size = 8, colour = "red", fontface = 2) +
  annotate("text", label = "DVI: 1.05", x = 1, y = -3.05, size = 4, colour = "black") +
  annotate("text", label = "DVI: 1.05", x = 2, y = -3.05, size = 4, colour = "black") +
  annotate("text", label = "DVI: 1.05", x = 3, y = -3.05, size = 4, colour = "black") +
  annotate("text", label = "DVI: 1.05", x = 4, y = -3.05, size = 4, colour = "black") +
  annotate("text", label = "DVI: 1.5", x = 5, y = -3.05, size = 4, colour = "black") +
  annotate("text", label = "DVI: 3.6", x = 6, y = -3.05, size = 4, colour = "black") 

rects <- tibble(xstart = c(0.5, 2.5), xend = c(2.5, 6.5), col = c("red", "blue"))

TableGrobGCaMP2 <- TableGCaMP %>% 
  select(Gene, Group1, Group2, HolmCorrection) %>% 
  mutate(Comparison = paste(Group1, " vs ", Group2)) %>% 
  select(Gene, Comparison, HolmCorrection) %>%
  rename(Holm_Corrected_p_value = HolmCorrection) %>%
  mutate(Holm_Corrected_p_value = round(Holm_Corrected_p_value, 7)) %>%
  pivot_wider(names_from = Gene, values_from = Holm_Corrected_p_value) %>%
  mutate(Comparison = case_when(Comparison == "OR5AN1 Line 15  vs  OR1A1 V2.2" ~ "OR5AN1 v OR1A1",
                                Comparison == "WT  vs  OR1A1 V2.2" ~ "WT v OR1A1",
                                Comparison == "WT  vs  OR5AN1 Line 15" ~ "WT v OR5AN1")) %>%
  select(Comparison, Olfr596_603, Olfr690, Olfr510, Olfr1154, Olfr358, Olfr390) %>%
  arrange(desc(Comparison))

ToSubsetx <- compare_means(DDCq ~ Genotype, Sample_Mutation_NO_Gene_GCaMP_DDCt, method = "t.test") %>% 
  select(group1, group2, p) %>%
  arrange(p) %>%
  mutate(row = row_number(),
         HolmAdjustment = 0.05/(3 - row + 1), 
         PassHA = p < HolmAdjustment, 
         HolmCorrection = p * (3 - row + 1),
         PassHC = HolmCorrection < 0.05, 
         Pass_DunnSidakCorrection = p < 1 - (1 - 0.05)^(1/3)) %>% 
  select(group1, group2, p, HolmCorrection, PassHC, Pass_DunnSidakCorrection) %>% 
  rename(p_value = p, Pass_HolmCorrection = PassHC, Group1 = group1, Group2 = group2) %>%
  rename(Holm_Corrected_p_value = HolmCorrection) %>%
  mutate(Holm_Corrected_p_value = round(Holm_Corrected_p_value, 7)) %>%
  rename(Combined_Olfr = Holm_Corrected_p_value)%>%
  mutate(Comparison = paste(Group1, " vs ", Group2)) %>% 
  mutate(Comparison = case_when(Comparison == "OR5AN1 Line 15  vs  OR1A1 V2.2" ~ "OR5AN1 v OR1A1",
                                Comparison == "WT  vs  OR1A1 V2.2" ~ "WT v OR1A1",
                                Comparison == "WT  vs  OR5AN1 Line 15" ~ "WT v OR5AN1")) %>%
  arrange(desc(Comparison))

TableGrobGCaMP2$Combined_Olfrs <- ToSubsetx$Combined_Olfr
TableGrobGCaMP2 <- TableGrobGCaMP2 %>%
  tableGrob(rows=NULL)

TextGrobGCaMP2 <- text_grob(paste("Eighteen independent t.tests were performed. With alpha = 0.05, Holm corrections for multiple comparisons was calculated. ", 
                                  "The table displays the corrected p-values. OR5AN1 is significantly different from WT for 4/6 Olfr mRNA values.", 
                                  "OR1A1 is significantly different from WT for 3/6 Olfr mRNA values. The non-WT lines differ statistically for Olfr690.",
                                  "Generating DDCt values by combining Olfr DCt together produces DDCt values that are significant versus WT, but not between 5AN1 and 1A1.", sep = "\n"), 
                            face = "italic", size = 10)

g <- gtable_add_grob(TableGrobGCaMP2,
                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                     t = 2, b = nrow(TableGrobGCaMP2), l = 1, r = ncol(TableGrobGCaMP2))

TableGrobGCaMP2 <- gtable_add_grob(g,
                                   grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                   t = 1, l = 1, r = ncol(g))

grid.arrange(L2FCGCaMP_2, TableGrobGCaMP2, TextGrobGCaMP2, heights = c(1.0, 0.18, 0.12), as.table = TRUE)

grid.arrange(L2FCGCaMP_3, TableGrobGCaMP2, heights = c(1.0, 0.18), as.table = TRUE)

#This one for that combines all the Olfr Genes into a single pool and calculates relevant stats and errors
Sample_Mutation_NO_Gene_GCaMP_DDCt <- Final_GCaMP_V2_Sample_MeanWT_DDCt %>%
  left_join(Final_GCaMP_V2_Mutation_WTmean_NoGeneDDCt, by = c("Mutation2")) %>%
  filter(Gene %notin% c("Acsm4", "Acss2", "Slc25a35")) %>%
  rename(Genotype = Mutation2) %>%
  mutate(Gene = factor(Gene, levels = c("Acsm4", "Acss2", "Slc25a35", "Olfr358", "Olfr390", "Olfr510", "Olfr596_603", "Olfr609", "Olfr690", "Olfr1154")),
         Genotype = case_when(Genotype == "L15_5AN1" ~ "OR5AN1 Line 15",
                              Genotype == "OR1A1_2.2" ~ "OR1A1 V2.2",
                              TRUE ~ "WT"))

Sample_Mutation_NO_Gene_GCaMP_DDCt %>%
  mutate(Genotype = case_when(Genotype == "OR5AN1 Line 15" ~ "OR5AN1", 
                              Genotype == "OR1A1 V2.2" ~ "OR1A1", 
                              TRUE ~ "WT")) %>%
  mutate(Genotype = factor(Genotype, levels = c("WT", "OR5AN1", "OR1A1"))) %>%
  ggplot(aes(Genotype, DDCq, color = Gene, group = Gene)) +
  geom_point(size = 3, alpha = 0.5) + 
  geom_errorbar(aes(ymin = DDCqMean - DDCqCI1prop, ymax = DDCqMean + DDCqCI1prop), width = 0.2) +
  labs(title = "Log2 Fold Change in Ct values for Grouped Olfr mRNA in OR5AN1 and OR1A1", 
       subtitle = "Plotted L2FC of each Sample with the 95% confidence interval based on the mean and single-propagated error of each Genotype", 
       caption = "DDCt calculated using corrected primer efficiencies and reference primers Acsm4, Acss2, and Slc25a35. \nEach Sample DDCt calibrated using the Gene-specific mean WT DCt. Genotype DDCt calibrated using mean WT Dct. Ct and DCt in triplicate. \nDisplayed are the Holm-corrected p values for three T-tests. The Mutant lines vs WT pass with Dunn-Šidák corrected alpha value. NS- Not Significant.", 
       y = "Log2 Fold Change (DDCt)") +
  geom_signif(y_position =c(1.4, 1.0, 0.6), comparisons = list(c("WT", "OR5AN1"), c("WT", "OR1A1"), c("OR5AN1", "OR1A1")), annotation = c("1.28e-6", "9.43e-10", "NS"), color = "black") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(face = "bold", size = 10),
        legend.position = c(0.92, 0.79),
        legend.background = element_rect(fill = "white", color = "black"), 
        plot.caption = element_text(size = 7, face = "italic"), 
        plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5))


#Just for reference to above
ToSubsetx <- compare_means(DDCq ~ Genotype, Sample_Mutation_NO_Gene_GCaMP_DDCt, method = "t.test") %>% 
  select(group1, group2, p) %>% 
  mutate(row = row_number(),
         HolmAdjustment = 0.05/(3 - row + 1), 
         PassHA = p < HolmAdjustment, 
         HolmCorrection = p * (3 - row + 1),
         PassHC = HolmCorrection < 0.05, 
         Pass_DunnSidakCorrection = p < 1 - (1 - 0.05)^(1/3)) %>% 
  select(group1, group2, p, HolmCorrection, PassHC, Pass_DunnSidakCorrection) %>% 
  rename(p_value = p, Pass_HolmCorrection = PassHC, Group1 = group1, Group2 = group2)
#----Exporting dataset----

#Cleaning up the Raw data file
Interm1 <- OR5AN1_6v15_Full %>%
  filter(Mutation %in% c("L15_WT", "L15_5AN1", "L6_WT")) %>%
  select(Plate, Mutation, Sample, Gene, Replicate, Cq) %>%
  mutate(Plate = case_when(Plate == "One" ~ "Seven",
                           Plate == "Two" ~ "Eight",
                           Plate == "Three" ~ "Nine"))

Raw_Data_1A1v5AN1 <- Final_OR1A1 %>%
  filter(Mutation %in% c("WT_2.2", "WT_1.1.2", "OR1A1_2.2")) %>%
  filter(Gene %notin% "Olfr609") %>%
  select(Plate, Mutation, Sample, Gene, Replicate, Cq) %>%
  mutate(Replicate = case_when(Replicate == 1 ~ "One",
                               Replicate == 2 ~ "Two",
                               Replicate == 3 ~ "Three")) %>%
  bind_rows(Interm1)

Raw_Data_1A1v5AN1_1 <- Raw_Data_1A1v5AN1 %>%
  select(Sample, Mutation, Gene, Plate, Replicate, Cq) %>%
  mutate(Gene = factor(Gene, levels = c("Acsm4", "Acss2", "Slc25a35", "Olfr358", "Olfr390", "Olfr510", "Olfr596_603", "Olfr609", "Olfr690", "Olfr1154")),
         Plate = factor(Plate, levels = c("One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine")),
         Replicate = factor(Replicate, levels = c("One", "Two", "Three"))) %>%
  rename(Genotype = Mutation, 
         Ct = Cq) %>%
  mutate(OR = case_when(Genotype %in% c("WT_2.2", "WT_1.1.2", "L15_WT", "L6_WT") ~ "WT",
                        Genotype == "OR1A1_2.2" ~ "OR1A1",
                        TRUE ~ "5AN1")) %>%
  mutate(Mouse_Line = case_when(Genotype %in% c("L15_WT", "L15_5AN1") ~ "5AN1-V2.2 Line 15",
                                Genotype == "L6_WT" ~ "5AN1-V2.2 Line 6",
                                Genotype == "WT_1.1.2" ~ "1A1-V1.1.2",
                                TRUE ~ "1A1-V2.2")) %>%
  mutate(Sample = gsub("\\_.*", "", Sample)) %>%
  select(Sample, Mouse_Line, OR, Gene, Plate, Replicate, Ct) %>%
  rename(Genotype = OR)

openxlsx::write.xlsx(Raw_Data_1A1v5AN1_1, file = "RawData_1A1v5AN1.xlsx")

#The first data compaction step, taking the mean of triplicates on each plate
Intermediate_Data_1A1v5AN1_1 <- Raw_Data_1A1v5AN1_1 %>%
  group_by(Sample, Mouse_Line, Genotype, Gene, Plate) %>%
  summarize(MeanSampleCt = mean(Ct),
            SampleSD = sd(Ct))

openxlsx::write.xlsx(Intermediate_Data_1A1v5AN1_1, file = "IntermediateData1_1A1v5AN1.xlsx")

#Calculating the Efficiency-corrected DCt values for each sample/gene from each plate
Interruption1 <- Intermediate_Data_1A1v5AN1_1 %>%
  select(-SampleSD) %>%
  left_join(AA10, by = "Gene") %>%
  ungroup() %>%
  mutate(MeanSCtxlog2Eff = MeanSampleCt * log2(Primer_Efficiency))

Interruption2 <- Interruption1 %>%
  filter(Gene %in% c("Acsm4", "Acss2", "Slc25a35")) %>%
  select(Sample, Gene, Plate, MeanSCtxlog2Eff) %>%
  pivot_wider(names_from = Gene, values_from = MeanSCtxlog2Eff) %>%
  rename(MeanSCtxlog2Eff_Acsm4 = Acsm4, 
         MeanSCtxlog2Eff_Acss2 = Acss2,
         MeanSCtxlog2Eff_Slc25a35 = Slc25a35)

Plate_Sample_DCt <- Interruption1 %>% 
  left_join(Interruption2, by = c("Sample", "Plate")) %>%
  mutate(DCtSample = (MeanSCtxlog2Eff_Acsm4 + MeanSCtxlog2Eff_Acss2 + MeanSCtxlog2Eff_Slc25a35)/3 - MeanSCtxlog2Eff)

openxlsx::write.xlsx(Plate_Sample_DCt, file = "Plate_Sample_DCt_1A1v5AN1.xlsx")

#Taking the mean of the DCtSample values from each plate to unify across plates
Sample_MeanDCt_SEM <- Plate_Sample_DCt %>%
  ungroup() %>% 
  group_by(Sample, Mouse_Line, Genotype, Gene) %>% 
  summarize(MeanDCtSample = mean(DCtSample), 
            DCtSampleSEM = sd(DCtSample)/sqrt(length(DCtSample)))

openxlsx::write.xlsx(Sample_MeanDCt_SEM, file = "Sample_MeanDCt_SEM_1A1v5AN1.xlsx")

#Calculating a DCt for each genotype/gene pair
Genotype_MeanDCt_SEM <- Sample_MeanDCt_SEM %>%
  ungroup() %>% 
  group_by(Genotype, Gene) %>% 
  summarize(MeanDCtGenotype = mean(MeanDCtSample),
            DCtGenotypeSize = length(MeanDCtSample),
            DCtGenotypeSEM = sd(MeanDCtSample)/sqrt(DCtGenotypeSize))

openxlsx::write.xlsx(Genotype_MeanDCt_SEM, file = "Genotype_MeanDCt_SEM_1A1v5AN1.xlsx")

#Calculating the DDCt for each genotype/gene pair
WTGenotype_MeanDCt_SEM <- Genotype_MeanDCt_SEM %>%
  filter(Genotype == "WT") %>%
  rename(DCtWTSize = DCtGenotypeSize, 
         MeanDCtWT = MeanDCtGenotype, 
         DCtWTSEM = DCtGenotypeSEM) %>%
  ungroup() %>%
  select(-Genotype)

DDCt_Genotype_SEM_CI <- Genotype_MeanDCt_SEM %>%
  left_join(WTGenotype_MeanDCt_SEM, by = "Gene") %>%
  mutate(MeanDDCt = MeanDCtGenotype - MeanDCtWT, 
         DDCtSEM = sqrt(DCtGenotypeSEM^2 + DCtWTSEM^2),
         DDCtCI = qt(0.975, df = DCtGenotypeSize + DCtWTSize - 2) * DDCtSEM)

openxlsx::write.xlsx(DDCt_Genotype_SEM_CI, file = "DDCt_Genotype_SEM_CI_1A1v5AN1.xlsx")

#DDCt for each sample; again for plotting purposes
DDCt_Sample_forPlots <- Sample_MeanDCt_SEM %>%
  ungroup() %>%
  left_join(WTGenotype_MeanDCt_SEM, by = "Gene") %>%
  select(Sample, Gene, Genotype, MeanDCtSample, MeanDCtWT) %>%
  mutate(SampleDDCt = MeanDCtSample - MeanDCtWT)

openxlsx::write.xlsx(DDCt_Sample_forPlots, file = "DDCt_Sample_forPlots_1A1v5AN1.xlsx")
#----More stats----
ForStats <- DDCt_Sample_forPlots %>%
  filter(Gene %notin% c("Acsm4", "Acss2", "Slc25a35")) #Removing the Reference Gene values from the table

#Let's explore this data set with the idea that I have three Genotypes and a column of DDCt values for each of them
#Let's plot the data points as a histogram for each Genotype
ggplot(ForStats, aes(SampleDDCt, fill = Gene)) + geom_histogram(binwidth = 0.5) +  facet_wrap(~Genotype)
ggplot(ForStats, aes(SampleDDCt, Genotype)) + geom_boxplot()

#Let's run a one-way ANOVA and then check the residuals
AOV1 <- aov(SampleDDCt ~ Genotype, data = ForStats) %>% summary()
#A histogram for the residuals
ggplot(ForStats, aes(AOV1$residuals)) + geom_histogram(binwidth = 0.25) + facet_wrap(~Genotype)
#A multiplot of the residuals and other metrics
par(mfrow=c(2,2))
plot(AOV1, lwd = 2)
par(mfrow=c(1,1))
#leveneTest and barlete test for homoscedasticity
library(car)
leveneTest(SampleDDCt ~ Genotype, data = ForStats)

bartlett.test(SampleDDCt ~ Genotype, data = ForStats)

fligner.test(SampleDDCt ~ Genotype, data = ForStats)

#tests for normality
shapiro.test(AOV1$residuals)

#non-parametric one-way anova 
kruskal.test(SampleDDCt ~ Genotype, data = ForStats)

#multiple comparisons using wilcox and t.test
compare_means(SampleDDCt ~ Genotype, data = ForStats, method = "t.test")
compare_means(SampleDDCt ~ Genotype, data = ForStats, method = "wilcox.test")

compare_means(SampleDDCt ~ Genotype, data = ForStats, method = "t.test") %>%
  arrange(p) %>%
  mutate(row = row_number(),
         HolmCorrection = p * (3 - row + 1),
         Pass_DunnSidakCorrection = p < 1 - (1 - 0.05)^(1/3))

#Now let's run Gene-split comparisons
ggplot(ForStats, aes(SampleDDCt)) + geom_histogram(binwidth = 0.25) +  facet_grid(rows = vars(Gene), cols = vars(Genotype))
ggplot(ForStats, aes(SampleDDCt, Genotype)) + geom_boxplot() + facet_wrap(~Gene)

#A quick ANOVA and krusal test for each gene group
aov(SampleDDCt ~ Genotype * Gene, data = ForStats) %>% summary()

compare_means(SampleDDCt ~ Genotype, data = ForStats, method = "anova", group.by = c("Gene"))
compare_means(SampleDDCt ~ Genotype, data = ForStats, method = "kruskal.test", group.by = c("Gene"))

#Making subgroups and running the toughest homoscedasticity test 
ForStats_358 <- ForStats %>% filter(Gene == "Olfr358")
ForStats_390 <- ForStats %>% filter(Gene == "Olfr390")
ForStats_510 <- ForStats %>% filter(Gene == "Olfr510")
ForStats_596_603 <- ForStats %>% filter(Gene == "Olfr596_603")
ForStats_690 <- ForStats %>% filter(Gene == "Olfr690")
ForStats_1154 <- ForStats %>% filter(Gene == "Olfr1154")

fligner.test(SampleDDCt ~ Genotype, data = ForStats_358)
fligner.test(SampleDDCt ~ Genotype, data = ForStats_390)
fligner.test(SampleDDCt ~ Genotype, data = ForStats_510)
fligner.test(SampleDDCt ~ Genotype, data = ForStats_596_603)
fligner.test(SampleDDCt ~ Genotype, data = ForStats_690)
fligner.test(SampleDDCt ~ Genotype, data = ForStats_1154)

bartlett.test(SampleDDCt ~ Genotype, data = ForStats_358)
bartlett.test(SampleDDCt ~ Genotype, data = ForStats_390)
bartlett.test(SampleDDCt ~ Genotype, data = ForStats_510)
bartlett.test(SampleDDCt ~ Genotype, data = ForStats_596_603)
bartlett.test(SampleDDCt ~ Genotype, data = ForStats_690)
bartlett.test(SampleDDCt ~ Genotype, data = ForStats_1154)

leveneTest(SampleDDCt ~ Genotype, data = ForStats_358)
leveneTest(SampleDDCt ~ Genotype, data = ForStats_390)
leveneTest(SampleDDCt ~ Genotype, data = ForStats_510)
leveneTest(SampleDDCt ~ Genotype, data = ForStats_596_603)
leveneTest(SampleDDCt ~ Genotype, data = ForStats_690)
leveneTest(SampleDDCt ~ Genotype, data = ForStats_1154)

#ANOVA for residuals for normality
AOV358 <- aov(SampleDDCt ~ Genotype, data = ForStats_358)
shapiro.test(AOV358$residuals)
AOV390 <- aov(SampleDDCt ~ Genotype, data = ForStats_390)
shapiro.test(AOV390$residuals)
AOV510 <- aov(SampleDDCt ~ Genotype, data = ForStats_510)
shapiro.test(AOV510$residuals)
AOV596_603 <- aov(SampleDDCt ~ Genotype, data = ForStats_596_603)
shapiro.test(AOV596_603$residuals)
AOV690 <- aov(SampleDDCt ~ Genotype, data = ForStats_690)
shapiro.test(AOV690$residuals)
AOV1154 <- aov(SampleDDCt ~ Genotype, data = ForStats_1154)
shapiro.test(AOV1154$residuals)

#Performing parametric and non-parametric pairwise comparisons
compare_means(SampleDDCt ~ Genotype, data = ForStats, method = "wilcox.test", group.by = c("Gene"))

compare_means(SampleDDCt ~ Genotype, data = ForStats, method = "t.test", group.by = c("Gene")) %>%
  select(Gene, group1, group2, p) %>%
  arrange(p) %>%
  mutate(row = row_number(),
         HolmAdjustment = 0.05/(18 - row + 1), 
         PassHA = p < HolmAdjustment, 
         HolmCorrection = p * (18 - row + 1),
         PassHC = HolmCorrection < 0.05, 
         Pass_DunnSidakCorrection = p < 1 - (1 - 0.05)^(1/18)) %>%
  mutate(NEWp = p.adjust(p, method = "holm"))

#ExTra plots
par(mfrow=c(2,2))
plot(AOV358, lwd = 2)
plot(AOV390, lwd = 2)
plot(AOV510, lwd = 2)
plot(AOV596_603, lwd = 2)
plot(AOV690, lwd = 2)
plot(AOV1154, lwd = 2)
par(mfrow=c(1,1))

#----CV values for this data----
ForCV1 <- Final_OR1A1 %>%
  filter(Mutation %in% c("WT_2.2", "WT_1.1.2", "OR1A1_2.2")) %>%
  filter(Gene %notin% "Olfr609") %>%
  select(Plate, Gene, Cq, Sample, Mutation) 

ForCV2 <- OR5AN1_6v15_Full %>%
  filter(Mutation %in% c("L15_WT", "L15_5AN1", "L6_WT")) %>%
  select(Plate, Gene, Cq, Sample, Mutation) %>%
  mutate(Plate = case_when(Plate == "One" ~ "A",
                           Plate == "Two" ~ "B",
                           Plate == "Three" ~ "C"))

CV3 <- bind_rows(ForCV1, ForCV2)

##Filtered Raw value Intra-Assay CV: 0.636
CV3 %>% 
  ungroup() %>% 
  select(Plate, Sample, Gene, Cq) %>%
  distinct() %>%
  group_by(Plate, Sample, Gene) %>%
  summarize(Plate = Plate, 
            Sample = Sample, 
            Gene = Gene, 
            RepMean = mean(Cq, na.rm = TRUE),
            RepSd = sd(Cq, na.rm = TRUE)) %>% 
  ungroup() %>% 
  distinct() %>%
  select(Plate, Sample, Gene, RepMean, RepSd) %>% 
  distinct() %>%
  mutate(CV = RepSd/RepMean * 100) %>%
  group_by(Plate) %>%
  summarize(MeanCV = mean(CV, na.rm = TRUE)) %>%
  summarize(mean = mean(MeanCV))


##Filtered Raw value Inter-Assay CV: 0.525
CV3 %>% 
  ungroup() %>% 
  select(Plate, Sample, Gene, Cq) %>%
  distinct() %>%
  group_by(Plate, Sample, Gene) %>%
  summarize(Plate = Plate, 
            Sample = Sample, 
            Gene = Gene, 
            RepMean = mean(Cq, na.rm = TRUE),
            RepSd = sd(Cq, na.rm = TRUE)) %>% 
  ungroup() %>% 
  select(Plate, Sample, Gene, RepMean, RepSd) %>% 
  distinct() %>%
  group_by(Sample, Gene) %>%
  summarize(PlateMeans = mean(RepMean, na.rm = TRUE), 
            PlateSd = sd(RepMean, na.rm = TRUE)) %>% 
  ungroup() %>% 
  summarise(CV = PlateSd/PlateMeans * 100) %>%
  summarize(meanCV = mean(CV))
