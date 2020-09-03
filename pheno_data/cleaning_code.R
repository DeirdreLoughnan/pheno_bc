rm(list=ls()) 
options(stringsAsFactors = FALSE)

setwd("~/Documents/github/pheno_bc/pheno_data")

#d<-read.csv("bb_phenobc_obs_Jan25.csv")
d<-read.csv("test_cleaning.csv", na.strings = "")

d
write.csv(d,"testoutput.csv")
names(d)

d$label<-paste(d$population,d$treatment,d$flask,d$species, sep="_")

d$bbch.t[d$label == "mp_LHH_NA_1_amealn"] <- 10000

head(d)

lateral<-list("mp_HC_LP_LF_26_samrac","mp_HC_HP_LF_21_samrac","mp_HC_HP_LF_10_alninc","
sm_HC_HP_LF_3_betpap","
sm_HC_HP_LF_14_poptre","
mp_HC_HP_LF_25_acegla","
mp_HC_LP_LF_39_poptre","
mp_HC_LP_LF_39_samrac","
mp_HC_LP_LF_36_popbal","
mp_LHH_24_poptre","
sm_LHH_25_poptre","
mp_LHH_32_poptre","
sm_LHH_37_poptre","
sm_LHH_37_sorsco","
sm_LHH_39_sorsco","
mp_LHL_1_poptre","
mp_LHL_2_poptre","
mp_LHL_3_spibet","
mp_LHL_4_popbal","
mp_LHL_8_popbal","
sm_LHL_9_popbal","
mp_LHL_14_poptre","
mp_LHL_18_alninc","
mp_LHL_24_popbal","
mp_LHL_26_poptre","
mp_LHL_27_poptre","
sm_LHL_34_poptre","
mp_LHL_36_acegla","
mp_LHL_40_alninc","
sm_LLH_18_popbal","
sm_LLH_25_sorsco","
mp_LLH_34_poptre","
sm_LLH_36_poptre","
sm_LLL_2_popbal","
mp_LLL_4_sorsco","
mp_LLL_7_popbal","sm_LLL_9_sorsco","mp_LLL_10_popbal","sm_LLL_10_alninc","mp_LLL_14_popbal","sm_LLL_19_sorsco","sm_LLL_24_alninc","sm_LLL_26_betpap","sm_LLL_29_popbal","mp_LLL_33_popbal","mp_LLL_35_poptre","mp_LLL_36_samrac","mp_LLL_40_poptre")


#mp_HC_LP_LF_26_samrac --> NA lateral
# mp_HC_HP_LF_21_samrac --> NA lateral
# mp_HC_HP_LF_10_alninc --> NA lateral
# sm_HC_HP_LF_3_betpap --> NA lateral
#sm_HC_HP_LF_14_poptre --> NA lateral
# mp_HC_HP_LF_25_acegla --> NA lateral
# mp_HC_LP_LF_39_poptre --> NA lateral
# mp_HC_LP_LF_39_samrac --> NA lateral
# mp_HC_LP_LF_36_popbal --> NA lateral
# mp_LHH_24_poptre --> NA lateral
# sm_LHH_25_poptre --> NA lateral
# mp_LHH_32_poptre --> NA lateral
# sm_LHH_37_poptre --> NA lateral
# sm_LHH_37_sorsco --> NA lateral
# sm_LHH_39_sorsco --> NA lateral
# mp_LHL_1_poptre --> NA lateral
# mp_LHL_2_poptre --> NA lateral
# mp_LHL_3_spibet --> NA lateral
# mp_LHL_4_popbal --> NA lateral
# mp_LHL_8_popbal --> NA lateral
# sm_LHL_9_popbal --> NA lateral
# mp_LHL_14_poptre --> NA lateral
# mp_LHL_18_alninc --> NA lateral
# mp_LHL_24_popbal --> NA lateral
# mp_LHL_26_poptre --> NA lateral
# mp_LHL_27_poptre --> NA lateral
# sm_LHL_34_poptre --> NA lateral
# mp_LHL_36_acegla --> NA lateral
# mp_LHL_40_alninc --> NA lateral
# sm_LLH_18_popbal --> NA lateral
# sm_LLH_25_sorsco --> NA lateral
# mp_LLH_34_poptre --> NA lateral
# sm_LLH_36_poptre --> NA lateral
# sm_LLL_2_popbal --> NA lateral
# mp_LLL_4_sorsco --> NA lateral
# mp_LLL_7_popbal --> NA lateral
# sm_LLL_9_sorsco --> NA lateral
# mp_LLL_10_popbal --> NA lateral
# sm_LLL_10_alninc --> NA lateral
# mp_LLL_14_popbal --> NA lateral
# sm_LLL_19_sorsco --> NA lateral
# sm_LLL_24_alninc --> NA lateral
# sm_LLL_26_betpap --> NA lateral
# sm_LLL_29_popbal --> NA lateral
# mp_LLL_33_popbal --> NA lateral
# mp_LLL_35_poptre --> NA lateral
# mp_LLL_36_samrac --> NA lateral
# mp_LLL_40_poptre --> NA lateral
 
# mp_LHH_9_samrac --> NA terminal bud
# sm_LHH_14_alninc --> NA terminal bud
# sm_LHH_15_spipyr --> NA terminal bud
# sm_LHH_18_amealn --> NA terminal bud
# mp_LHH_25_alninc --> NA terminal bud
# mp_LHL_14_loninv --> NA terminal bud
# sm_LHL_14_rubpar --> NA terminal bud
# mp_LHL_15_rubpar --> NA terminal bud
# sm_LHL_19_riblac --> NA terminal
# mp_LHL_25_spipyr --> NA terminal bud
# sm_LLH_16_vacmem --> NA terminal bud
# sm_LLH_19_poptre --> NA terminal bud
# sm_LLH_21_popbal --> NA terminal bud
# sm_LLH_28_loninv --> NA terminal bud
# mp_LLL_3_shecan --> NA terminal bud
# sm_LLL_5_betpap --> NA terminal bud
# sm_LLL_5_popbal --> NA terminal bud
# sm_LLL_7_loninv --> NA terminal bud
# mp_LLL_30_sorsco --> NA terminal bud
# sm_LLL_32_alnvir --> NA terminal bud
# mp_LLL_22_corsto --> NA terminal bud
# mp_LLL_34_popbal --> NA terminal


#dead: 
#sm_LLH_23_amealn
#mp_LLL_28_acegla
#sm_HHH_10_rupar
#replace all symalb terminal buds with NA

#replace all menfer terminal with NA, what I thought water terminal buds were flowers