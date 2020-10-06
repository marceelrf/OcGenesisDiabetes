# OcGenesisDiabetes
Osteoclastogenesis genes in diabetes

### Data sets
#### GSE54154
+ BMGM_2: Diabetes do tipo 2 diferenciados com GM-CSF (BMGM);
+ BMGM_db: Não diabeticos diferenciados com M-CSF (BMMC);
+ BMMC_2: Diabetes do tipo 2 diferenciados com M-CSF (BMMC);
+ BMMC_db: Não diabeticos diferenciados com GM-CSF (BMGM);

GM-CSF usado para diferenciar em macrófago do tipo 1 (M1);  
M-CSF usado para diferenciar em macrófago do tipo 2 (M2);


##### Constrastes analisados
A análise de expressão diferencial foi realizada utilizando o pacote edgeR.  
Constrastes:  
BMMC_db_vs_BMGM_db = BMMC_db - BMGM_db #M2 vs M1 - Controle
BMMC_2_vs_BMGM_2 = BMMC_2 - BMGM_2 #M2 vs M1 - Diabetes tipo 2
BMGM_db_vs_BMGM_2 = BMGM_db - BMGM_2 #M1 controle vs M1 diabetes tipo 2
BMMC_db_vs_BMMC_2 = BMMC_db - BMMC_2 #M2 controle vs M2 diabetes tipo 2

Arquivo `DEGs` contem os genes diferencialmente expressos baseados nos contrastes.
#### E-GEOD-54779
