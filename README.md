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
