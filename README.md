# Antoine.code

```{r}
library(dada2); packageVersion("dada2") # importation de dada2
```



```{r}
path <- "~/Analyse Article ADM/MiSeq_SOP" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)  # Cette fonction permet de lister le fichier fastq
```
 [1] "F3D0_S188_L001_R1_001.fastq"                          
 [2] "F3D0_S188_L001_R2_001.fastq"                          
 [3] "F3D1_S189_L001_R1_001.fastq"                          
 [4] "F3D1_S189_L001_R2_001.fastq"                          
 [5] "F3D141_S207_L001_R1_001.fastq"                        
 [6] "F3D141_S207_L001_R2_001.fastq"                        
 [7] "F3D142_S208_L001_R1_001.fastq"                        
 [8] "F3D142_S208_L001_R2_001.fastq"                        
 [9] "F3D143_S209_L001_R1_001.fastq"                        
[10] "F3D143_S209_L001_R2_001.fastq"                        
[11] "F3D144_S210_L001_R1_001.fastq"                        
[12] "F3D144_S210_L001_R2_001.fastq"                        
[13] "F3D145_S211_L001_R1_001.fastq"                        
[14] "F3D145_S211_L001_R2_001.fastq"                        
[15] "F3D146_S212_L001_R1_001.fastq"
[16] "F3D146_S212_L001_R2_001.fastq"                        
[17] "F3D147_S213_L001_R1_001.fastq"                        
[18] "F3D147_S213_L001_R2_001.fastq"                        
[19] "F3D148_S214_L001_R1_001.fastq"                        
[20] "F3D148_S214_L001_R2_001.fastq"                        
[21] "F3D149_S215_L001_R1_001.fastq"                        
[22] "F3D149_S215_L001_R2_001.fastq"                        
[23] "F3D150_S216_L001_R1_001.fastq"                        
[24] "F3D150_S216_L001_R2_001.fastq"                        
[25] "F3D2_S190_L001_R1_001.fastq"                          
[26] "F3D2_S190_L001_R2_001.fastq"                          
[27] "F3D3_S191_L001_R1_001.fastq"                          
[28] "F3D3_S191_L001_R2_001.fastq"                          
[29] "F3D5_S193_L001_R1_001.fastq"                          
[30] "F3D5_S193_L001_R2_001.fastq"                          
[31] "F3D6_S194_L001_R1_001.fastq"                          
[32] "F3D6_S194_L001_R2_001.fastq" 
[33] "F3D7_S195_L001_R1_001.fastq"                          
[34] "F3D7_S195_L001_R2_001.fastq"                          
[35] "F3D8_S196_L001_R1_001.fastq"                          
[36] "F3D8_S196_L001_R2_001.fastq"                          
[37] "F3D9_S197_L001_R1_001.fastq"                          
[38] "F3D9_S197_L001_R2_001.fastq"                          
[39] "filtered"                                             
[40] "HMP_MOCK.v35.fasta"                                   
[41] "Mock_S280_L001_R1_001.fastq"                          
[42] "Mock_S280_L001_R2_001.fastq"                          
[43] "mouse.dpw.metadata"                                   
[44] "mouse.time.design"                                    
[45] "silva_nr99_v138.2_toGenus_trainset.fa.gz?download=1"  
[46] "silva_nr99_v138.2_toSpecies_trainset.fa.gz?download=1"
[47] "SILVA_SSU_r138_2_2024.RData"                          
[48] "silva_v138.2_assignSpecies.f
[49] "stability.batch"                                      
[50] "stability.files" 

On commence par indiquer le dossier où sont stockés les fichiers FASTQ, puis on affiche son contenu afin de vérifier la présence des échantillons et des fichiers de référence.

