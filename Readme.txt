# Read me
The data underlying this article are available at https://github.com/LingxuanShao/Robust-FDA.

(1) The files "estimation_rho1_tuned_kappa.R" and "estimation_rho2_rho3.R" contain the code used for assessing our method in Section 5.1.

(2) The file "estimation_score.R" contains the code used for score estimation and comparison with PACE in Section 5.2.

(3) The file "AZcode.R" contains the code used for the real data analysis in Section 6.

(4) The file "hippo.RDATA" provides the Alzheimer's data: "metainfo" is a data.frame of size 977*14, containing information for each observation. Specifically, metainfo[i,4], metainfo[i,5], and metainfo[i,7] represent the subject ID, the signal for CN/AD, and the observed time for the i-th observation, respectively. "dti" is an array that includes the FA data. For further details, visit https://github.com/linulysses/iRFDA-sparse.

