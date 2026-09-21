
#  TFM          25/09/2026 
#     
#  Ana Pilar Melendo Vidal
# ==============================================================================

#==========================
# Se cargan los archivos
#==========================

metilacion <- read.delim("CCLE_RRBS_TSS_1kb_20180614.txt", header=T, sep="\t", fill=T, row.names=1)
matriz <- read.csv("OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv", sep = ",", header = TRUE, row.names = NULL)
model <- read.csv("Model.csv", sep = ",", header = TRUE, row.names = NULL)
cn <- read.csv("OmicsCNGeneWGS.csv", sep = ",", header = TRUE, row.names = NULL)
mutacion <- read.csv("OmicsSomaticMutations.csv", sep = ",", header = TRUE, row.names = NULL)
farma1 <- read.csv("Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv", sep = ",", header = TRUE, row.names = NULL)
farma2 <- read.csv("Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv", sep = ",", header = TRUE, row.names = NULL)
MSI <- read.csv("OmicsGlobalSignatures.csv")

#===========================
# Se cargan las librerias
#===========================

library(org.Hs.eg.db)
library(AnnotationDbi)
library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(car)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(clusterProfiler)
library(glmnet)
library(foreach)    
library(doParallel)
library(doRNG) 


#==============================================================================
# FASE 1: ModelID comunes, taxonomía tisular e integración de cofactores
#==============================================================================

# -------------------------------------------
# 1.1. ModelID comunes de fármacos y ómicas
# -------------------------------------------

model <- model %>% filter(OncotreePrimaryDisease != "Non-Cancerous") #Se eliminan líneas de control

farma1 <- farma1 %>%
  rename(IDs = 1) %>% # Se cambia el nombre de la primera columna
  rename_with(~ gsub("\\.", "-", .x)) # Se cambian los puntos por guiones


ids_farma <- colnames(farma1)[startsWith(colnames(farma1), "ACH-")] # Se extraen los IDs de farma1

ids_expr <- unique(matriz$ModelID[matriz$IsDefaultEntryForModel == "Yes"]) # Se extraen los IDs de expresión
 
ccle_names_meth <- colnames(metilacion)[!colnames(metilacion) %in% c("gene", "chr", "fpos", "tpos", "strand", "avg_coverage")]
ids_meth  <- unique(model$ModelID[model$CCLEName %in% ccle_names_meth])# Se extraen los IDs de metilación

model_ids_comunes <- intersect(intersect(ids_farma, ids_expr), ids_meth)# Líneas celulares comunes iniciales

#-------------------------------------------
# 1.2. Filtrado y traducción de linajes
#-------------------------------------------

linajes_excluir <- c("Normal", "Fibroblast", "Hair", "", "Other") #Se excluyen linajes de control o tumorales no definidos

lineas_excluidas_linaje_control <- model %>% filter(ModelID %in% model_ids_comunes) %>% filter(OncotreeLineage %in% linajes_excluir)
length(nrow(lineas_excluidas_linaje_control)) # Líneas excluidas por ser control

lineas_exclu_repre <- model %>% 
  filter(ModelID %in% model_ids_comunes) %>%
  filter(!OncotreeLineage %in% linajes_excluir) %>%  group_by(OncotreeLineage) %>%
  filter(n() < 15) %>% # Filtro de representatividad estadística
  ungroup() 


MSI <- MSI %>% filter(IsDefaultEntryForModel == "Yes") %>% dplyr::select(-X)


model_global <- model %>%
  filter(ModelID %in% model_ids_comunes) %>%
  filter(!OncotreeLineage %in% linajes_excluir) %>% 
  left_join(MSI, by = "ModelID") %>%  
  mutate(
    growth_binario = case_when(
      GrowthPattern == "Adherent"   ~ 1,
      GrowthPattern == "Mixed"      ~ 1, #Se incluyen con las adherentes
      GrowthPattern == "Suspension" ~ 0,
      TRUE                          ~ NA_real_
    ),
    model_type_factor = as.factor(ModelType), # Se factoriza el tipo de modelo celular
    # Traducción
    OncotreeLineage = case_when(
      OncotreeLineage == "Lung" ~ "Pulmón",
      OncotreeLineage == "Lymphoid" ~ "Linfoide",
      OncotreeLineage == "Esophagus/Stomach" ~ "Esófago/Estómago",
      OncotreeLineage == "Skin" ~ "Piel",
      OncotreeLineage == "CNS/Brain" ~ "SNC/Cerebro",
      OncotreeLineage == "Pancreas" ~ "Páncreas",
      OncotreeLineage == "Myeloid" ~ "Mieloide",
      OncotreeLineage == "Bowel" ~ "Intestino",
      OncotreeLineage == "Breast" ~ "Mama",
      OncotreeLineage == "Head and Neck" ~ "Cabeza y Cuello",
      OncotreeLineage == "Kidney" ~ "Riñón",
      OncotreeLineage == "Liver" ~ "Hígado",
      OncotreeLineage == "Bladder/Urinary Tract" ~ "Vejiga/Tracto Urinario",
      OncotreeLineage == "Bone" ~ "Hueso",
      OncotreeLineage == "Peripheral Nervous System" ~ "Sistema Nervioso Periférico",
      OncotreeLineage == "Ovary/Fallopian Tube" ~ "Ovario/Trompas de Falopio",
      OncotreeLineage == "Prostate" ~ "Próstata",
      OncotreeLineage == "Uterus" ~ "Útero",
      OncotreeLineage == "Cervix" ~ "Cérvix",
      OncotreeLineage == "Vulva/Vagina" ~ "Vulva/Vagina",
      OncotreeLineage == "Testis" ~ "Testículo",
      TRUE ~ as.character(OncotreeLineage)
    )) %>%
  filter(!is.na(OncotreeLineage) & !is.na(MSIScore) & !is.na(model_type_factor) & !is.na(growth_binario)) %>% 
  group_by(OncotreeLineage) %>%
  filter(n() >= 15) %>% # Filtro de representatividad estadística
  ungroup() %>% 
  dplyr::select(CCLEName, ModelID, OncotreeLineage, OncotreePrimaryDisease, growth_binario, model_type_factor, MSIScore) 

unique(model_global$OncotreeLineage) #Linajes que pasan el filtro


#=======================================================
# FASE 2: Modulo de epigenética (rango dinámico y PCA)
#=======================================================

# 2.1. Se eliminan los genes de los cromosomas X e Y
#-----------------------------------------------------

ccle_names_definitivos <- model_global$CCLEName

metilacion_sin_gsex <- metilacion [!metilacion$chr %in% c("chrX", "chrY"), ]

nfilas_meth_con_gsex <- length(rownames(metilacion)) # Genes iniciales

nfilas_meth_sin_gsex <- length(rownames(metilacion_sin_gsex)) # Genes tras eliminar los de los cromosomas X e Y


# 2.2. Preparación de la matriz numérica
#----------------------------------------

mat_meth_numeric <- metilacion_sin_gsex %>% 
  dplyr::select(all_of(intersect(ccle_names_meth, ccle_names_definitivos))) %>% 
  mutate(across(everything(), as.numeric)) %>% 
  as.matrix()

n_lineas_meth <- ncol(mat_meth_numeric) # Número de lineas celulares iniciales


#2.3. Selección de genes con datos de metilación en el 80% de las líneas celulares mínimo (20% de NAs máximo)
#-------------------------------------------------------------------------------------------------------------

prop_nas_por_cpg <- rowMeans(is.na(mat_meth_numeric))# Se calcula la proporción de NAs por fila (cada fila es un gen/sitio CpG)

n_cpgs_iniciales <- length(prop_nas_por_cpg) # Registro los CpGs de partida antes del filtro 


print(summary(prop_nas_por_cpg)) # Reporte de datos faltantes


umbral_nas <- 0.20 #Umbral del 20%

cpgs_conservados <- prop_nas_por_cpg <= umbral_nas

n_cpgs_filtrados_nas <- sum(!cpgs_conservados) # Promotores que se eliminan

mat_meth_numeric <- mat_meth_numeric[cpgs_conservados, , drop = FALSE] #Se eliminan los genes que no cumplen el umbral


#2.4. Sustitución de NAs por la mediana de ese gen/CpG
#-----------------------------------------------------

row_medians <- apply(mat_meth_numeric, 1, median, na.rm = TRUE) # Se calcula la mediana de metilación de cada fila (cada sitio CpG/gen en las diferentes lineas celulares)

na_indices <- which(is.na(mat_meth_numeric), arr.ind = TRUE)  # Se localizan los NAs de la matriz

mat_meth_numeric[na_indices] <- row_medians[na_indices[, 1]] # Se reemplaza cada NA por la mediana correspondiente a su fila (su sitio CpG)


mat_pca_meth <- t(mat_meth_numeric) # Se transpone la matriz (Filas = Muestras, Columnas = Genes/Sitios)

n_genes_partida_pca <- ncol(mat_pca_meth) # Total de genes (promotores) listos para análisis de varianza


# 2.5. PCA global
#-----------------

varianzas_finales <- apply(mat_pca_meth, 2, var, na.rm = TRUE) # Se calcula la varianza por columna

n_top_pca <- 15000 # Se seleccionan las 15.000 variables con mayor variabilidad biológica

top_15000_promotores <- names(sort(varianzas_finales, decreasing = TRUE)[1:min(n_top_pca, length(varianzas_finales))]) 

mat_pca_top15000 <- mat_pca_meth[, top_15000_promotores, drop = FALSE] 

pca_meth_res <- prcomp(mat_pca_top15000, scale. = TRUE) # Se calcula el PCA 

# Extracción de coordenadas
coordenadas_pca <- data.frame(
  CCLEName = rownames(pca_meth_res$x),
  PC1 = as.numeric(pca_meth_res$x[, 1]),
  PC2 = as.numeric(pca_meth_res$x[, 2])
) %>%
  inner_join(dplyr::select(model_global, CCLEName, ModelID, OncotreeLineage), by = "CCLEName") %>%
  distinct(ModelID, .keep_all = TRUE)


colores <- c(
  "Linfoide" = "#00FF00", "Mieloide" = "#FF00FF", 
  "Pulmón" = "#1F78B4", "Mama" = "#4A148C", "Intestino" = "#78909C", 
  "Esófago/Estómago" = "#A6CEE3", "SNC/Cerebro" = "#B0BEC5", "Páncreas" = "#CFD8DC", 
  "Piel" = "#8D6E63", "Cabeza y Cuello" = "#BCAAA4", "Riñón" = "#9FA8DA", 
  "Hígado" = "#5C6BC0", "Vejiga/Tracto Urinario" = "#81C784", "Hueso" = "#7570B3", 
  "Sistema Nervioso Periférico" = "#4DB6AC", "Ovario/Trompas de Falopio" = "#A1887F", 
  "Próstata" = "#D7CCC8", "Útero" = "#FFE0B2", "Cérvix" = "#B39DDB", 
  "Vulva/Vagina" = "#D1C4E9", "Testículo" = "#EEEEEE"
)

ggplot(coordenadas_pca, aes(x = PC1, y = PC2, color = OncotreeLineage)) +
  geom_point(alpha = 0.8, size = 2.5) +  
  scale_color_manual(values = colores) + 
  theme_classic(base_size = 22) +
  labs(
    x = "PC1",
    y = "PC2",
    color = "Tejido de origen"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 24),
    axis.title = element_text(face = "bold", size = 22),
    legend.title = element_text(face = "bold", size = 22),
    legend.text = element_text(size = 22))


matriz_pcs_limpia <- coordenadas_pca %>% dplyr::select(ModelID, PC1, PC2)

model_global <- model_global %>%
  inner_join(matriz_pcs_limpia, by = "ModelID") %>%
  dplyr::select(CCLEName, ModelID, OncotreeLineage, OncotreePrimaryDisease, GrowthPattern = growth_binario, MSIScore, PC1, PC2) %>%
  distinct(ModelID, .keep_all = TRUE) # Se quitan duplicados


# 2.6. Cribado por el rango dinámico epigenético (Delta P90 - P10 >= 0.4)
#-------------------------------------------------------------------------

rangos_dinamicos <- apply(mat_pca_meth, 2, function(col) {
  quantile(col, 0.90, na.rm = TRUE) - quantile(col, 0.10, na.rm = TRUE)})

genes_rango_dinamico <- names(rangos_dinamicos[rangos_dinamicos >= 0.40]) #Genes que cumplen el rango

mat_meth_supervivientes <- mat_pca_meth[, genes_rango_dinamico, drop = FALSE]

methyl_final <- as.data.frame(t(mat_meth_supervivientes)) %>%
  mutate(Gene = rownames(.)) %>%
  pivot_longer(
    cols = -Gene,
    names_to = "CCLEName",
    values_to = "metilacion"
  ) %>%
  inner_join(model_global, by = "CCLEName")

methyl_final$Gene <- as.character(methyl_final$Gene)

methyl_final$Gene <- sub("_.*", "", methyl_final$Gene)#Se modifica el nombre de los genes 


# 2.7. Reporte 
#--------------

cat("Líneas celulares analizadas: ", n_lineas_meth, "\n")
cat("Promotores (CpGs) iniciales evaluados: ", nfilas_meth_con_gsex, "\n")
cat("Pomotores resultantes al quitar genes de cromosomas X e Y: ", nfilas_meth_sin_gsex, "\n")
cat("Promotores descartados por exceso de NAs (>20%): ", n_cpgs_filtrados_nas, "\n")
cat("Promotores finales con rango dinámico: ", length(genes_rango_dinamico), "\n") 
cat("Eventos CpG-línea totales de partida:  ", (n_cpgs_iniciales * n_lineas_meth), "\n")
cat("Eventos CpG-línea finales: ", nrow(methyl_final), "\n")


#==================================================================
# FASE 3: Módulo del transcriptoma (umbral mínimo y de actividad)
#==================================================================

# 3.1. Preparación de la matriz
#--------------------------------
mat_clean <- matriz %>% 
  filter(IsDefaultEntryForModel == "Yes") %>% 
  filter(ModelID %in% model_global$ModelID) %>% #ModelsID comunes
  dplyr::select(ModelID, where(is.numeric)) %>%
  rename_with(~ str_replace(., "\\.\\..*", ""), -ModelID) %>%
  dplyr::select(-X)
  
n_lineas_finales  <- nrow(mat_clean) 

n_genes_iniciales <- ncol(mat_clean) - 1 


# 3.2. Cálculo del umbral mínimo de la densidad
#-----------------------------------------------
expr_long_previa <- mat_clean %>%
  pivot_longer(cols = -ModelID, names_to = "Gene", values_to = "TPM")

densidad_datos <- density(expr_long_previa$TPM, na.rm = TRUE, from = 0, to = 6)

tabla_densidad <- data.frame(x = densidad_datos$x, y = densidad_datos$y)

umbral_matematico_exacto <- tabla_densidad %>%
  filter(x > 0.5 & x < 2.5) %>% # Se acota la búsqueda
  filter(y == min(y)) %>% # Mínimo local
  pull(x)


ggplot(expr_long_previa, aes(x = TPM)) +
  geom_density(fill = "gray85", color = "gray30", linewidth = 0.5, alpha = 0.7) + 
  geom_vline(xintercept = umbral_matematico_exacto, color = "darkred", linetype = "dashed", linewidth = 1) + 
  labs(
    x = "Nivel de expresión [log2 (TPM+1)]", 
    y = "Densidad",
    title = "Distribución global del transcriptoma"
  ) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA), breaks = seq(0, 1, by = 0.25)) +  
  scale_x_continuous(expand = c(0, 0), limits = c(0, 15), breaks = seq(0, 14, by = 2)) + 
  coord_cartesian(xlim = c(0, 15), ylim = c(0, NA), clip = "off") +
  theme_classic(base_size = 22) + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 24), 
    axis.title = element_text(face = "bold", size = 22),
    axis.line = element_line(color = "black", linewidth = 0.5), 
    axis.ticks = element_line(color = "black")                  
  )


# 3.3. Filtro actividad (>= 10%)
# -------------------------------
mat_expr_numerica <- mat_clean %>% dplyr::select(-ModelID) %>% as.matrix()

proporciones_genes <- colSums(mat_expr_numerica > umbral_matematico_exacto, na.rm = TRUE) / n_lineas_finales #Porcentaje de líneas en las que el gene esta activo

genes_pasan_filtro <- names(proporciones_genes)[proporciones_genes >= 0.10]

expr_filtrada <- mat_clean %>%
  dplyr::select(ModelID, all_of(genes_pasan_filtro)) %>% #Se seleccionan solo los genes activos
  pivot_longer(cols = -ModelID, names_to = "Gene", values_to = "TPM")


# 3.4. Reporte 
#--------------

n_genes_descartados  <- n_genes_iniciales - (length(genes_pasan_filtro))

porcentaje_descarte  <- (n_genes_descartados / n_genes_iniciales) * 100

cat("Líneas celulares analizadas: ", n_lineas_finales, "\n")
cat("Umbral mínimo de la densidad calculado: TPM =", round(umbral_matematico_exacto, 3), "\n") 
cat("Genes iniciales evaluados: ", n_genes_iniciales, "\n")
cat("Genes activos (presencia >= 10%): ", length(genes_pasan_filtro), "\n")
cat("Genes filtrados por inactividad constitutiva: ", n_genes_descartados, "\n") 
cat("Tasa porcentual de descarte de genes: ", round(porcentaje_descarte, 2), "%\n")
cat("Eventos gen-línea totales de partida: ", (n_genes_iniciales * n_lineas_finales), "\n") 
cat("Eventos gen-línea activos finales (listos):  ", nrow(expr_filtrada), "\n")


#====================================================================================================
# FASE 4: Módulo de genoma (mutaciones somáticas LoF y deleciones)
#====================================================================================================

# 4.1. Procesamiento mutaciones somáticas (filtro LoF / NMD)
#------------------------------------------------------------

mutacion_limpia <- mutacion %>% 
  dplyr::rename(Gene = "HugoSymbol") %>% 
  filter(IsDefaultEntryForModel == "Yes") %>% 
  filter(LikelyLoF == "True" | grepl("NMD", VariantInfo)) %>% #Se eliminan genes que pueden reducir la expresión génica
  filter(ModelID %in% model_global$ModelID) %>% # ModelID comunes
  filter(Gene %in% genes_pasan_filtro) %>% #Genes activos    
  mutate(mut = TRUE) %>% 
  distinct(ModelID, Gene, mut)

n_mut_lineas_afectadas <- length(unique(mutacion_limpia$ModelID)) #Líneas celulares afectadas

n_genes_mutados_lof <- length(unique(mutacion_limpia$Gene)) #Genes mutados


# 4.2. Procesamiento del número de copias (CN <= 1)
#---------------------------------------------------------------
n_genes_cn_iniciales <- ncol(cn) - 6 # Conteo inicial de genes 

cn_limpio <- cn %>% 
  dplyr::select(-X) %>% 
  filter(IsDefaultEntryForModel == "Yes") %>% 
  filter(ModelID %in% model_global$ModelID) %>% # ModelID comunes
  dplyr::select(ModelID, where(is.numeric)) %>%   
  rename_with(~ str_replace(., "\\.\\..*", ""), -ModelID) %>%
  pivot_longer(
    cols = -ModelID,
    names_to = "Gene", 
    values_to = "CN"
  ) %>%
  filter(Gene %in% genes_pasan_filtro) %>% 
  distinct(ModelID, Gene, CN) 

cn_limpio1 <- cn_limpio %>% filter(CN <= 1) 
nrow(cn_limpio1) #Registros con deleciones


# 4.3. Aislamiento evento moleculares puros
#--------------------------------------------

eventos_transcriptomicos_entrada <- nrow(expr_filtrada) # Conteo eventos transcriptómicos antes del filtro

expr_mut_cn <- expr_filtrada %>% 
  left_join(mutacion_limpia, by = c("ModelID", "Gene")) %>% # Se une expresión con mutaciones
  left_join(cn_limpio, by = c("ModelID", "Gene")) %>% # Se une expresión con CN
  mutate(
    mut = if_else(is.na(mut), FALSE, mut), 
    del = if_else(is.na(CN), FALSE, CN <= 1) 
  ) %>%
  dplyr::select(-CN) %>% 
  filter(mut == FALSE & del == FALSE) %>%  # Se conservan genes activos sin mutaciones LoF ni CN
  dplyr::select(-mut, -del) 


# 4.4. Reporte
#--------------
eventos_rotos_eliminados <- eventos_transcriptomicos_entrada - (nrow(expr_mut_cn))
porcentaje_perdida_daño  <- (eventos_rotos_eliminados / eventos_transcriptomicos_entrada) * 100

dif1 <- as.data.frame (cn_limpio1)[, c("ModelID", "Gene")]
dif2 <- as.data.frame (mutacion_limpia)[, c("ModelID", "Gene")]
comuness <- intersect (dif1, dif2)

cat("Genes únicos con mutaciones LoF: ", n_genes_mutados_lof, "\n")
cat("Líneas celulares afectadas por eventos LoF/NMD: ", n_mut_lineas_afectadas, "\n") 
cat("Registros de expresión de entrada:", eventos_transcriptomicos_entrada, "\n") 
cat("Registros transcriptómicos puros:  ", nrow(expr_mut_cn), "\n") 
cat("Registros descartados por daño genómico (Mutación/Del): ", eventos_rotos_eliminados, "\n") 
cat("Tasa de pérdida de eventos por daño estructural: ", round(porcentaje_perdida_daño, 2), "%\n") 
cat("De los registros descartados, son solo mutaciones", nrow(mutacion_limpia), "\n")
cat("De los registros descartados, son solo deleciones", nrow(cn_limpio1), "\n")
cat("De los registros descartados, están duplicados (son mutaciones y deleciones)", nrow(comuness), "\n")


#========================================
# FASE 5: Evaluación de las covariables       
#========================================

# 5.1. Cruce multiómico
#------------------------

tabla_modelos <- expr_mut_cn %>%
  inner_join(methyl_final %>% dplyr::select(ModelID, CCLEName, Gene, metilacion), by = c("ModelID", "Gene")) %>%
  inner_join(model_global, by = c("ModelID", "CCLEName")) #Union eventos transcriptómicos puros con metadatos globales


# 5.2. Comprobación representatividad de linajes (n >= 15)
#---------------------------------------------------------

linajes_validos_finales <- tabla_modelos %>%
  distinct(ModelID, OncotreeLineage) %>%
  count(OncotreeLineage) %>%
  filter(n >= 15) %>%
  pull(OncotreeLineage) 

linajes_eliminados <- tabla_modelos %>%
  distinct(ModelID, OncotreeLineage) %>%
  count(OncotreeLineage) %>%
  filter(n < 15) %>%
  pull(OncotreeLineage)

tabla_modelos <- tabla_modelos %>% filter(OncotreeLineage %in% linajes_validos_finales)


# 5.3. Evaluación multicolinealidad (VIF)
#-----------------------------------------

datos_vif_real <- tabla_modelos %>%
  distinct(ModelID, .keep_all = TRUE) %>%
  mutate(GrowthPattern = as.factor(GrowthPattern))

modelo_piloto_vif_real <- lm(TPM ~ metilacion + OncotreeLineage + GrowthPattern + MSIScore + PC1 + PC2, data = datos_vif_real) # Modelo lineal piloto

tabla_vif_real <- car::vif(modelo_piloto_vif_real) # Cálculo de valores GVIF

print(tabla_vif_real) #Valores GVIF de las covariables


#===================================================
# FASE 6: Se define la población celular a analizar
#===================================================

lineas_totales_tfm <- length(unique(tabla_modelos$ModelID)) # Lineas celulares finales a analizar

genes_totales_tfm <- length(unique(tabla_modelos$Gene)) 

eventos_totales_tfm <- nrow(tabla_modelos)

prop.table(table(tabla_modelos$GrowthPattern)) * 100 # Distribución del patrón de crecimiento

summary(tabla_modelos$MSIScore)

metadatos_modelos <- model %>% inner_join(model_global, by = "ModelID")

prop.table(table(metadatos_modelos$PatientRace)) * 100 # Distribución de la etnia

prop.table(table(metadatos_modelos$Sex)) * 100 # Distribución del género
  
prop.table(table(metadatos_modelos$AgeCategory)) * 100 # Distribución de la edad
  
# Frecuencias de líneas celulares por tipo de cáncer
df_linajes_sincronizados <- tabla_modelos %>% 
  distinct(ModelID, OncotreeLineage) %>% 
  count(OncotreeLineage) %>% arrange(desc(n))

ggplot(df_linajes_sincronizados, aes(x = reorder(OncotreeLineage, n), y = n)) +
  geom_bar(stat = "identity", fill = "#4A5568", alpha = 0.8) + 
  geom_text(aes(label = n), hjust = -0.3, size = 5, fontface = "bold", color = "#1A1A1A") +
  coord_flip() + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(df_linajes_sincronizados$n) * 1.10)) + 
  theme_classic(base_size = 22) + 
  labs(
    title = "Líneas celulares por linaje tumoral",
    x = "Tejido de origen",
    y = "Número de muestras (líneas celulares)"
  ) +
  theme(
    plot.title = element_text(face = "bold", size=24),
    axis.title = element_text(face = "bold", size=22))


#====================================================================================
# FASE 7: Genes que reducen su expresión al aumentar la metilacion de sus promotores
#====================================================================================

# 7.1. Modelo mixto multivariante 
#----------------------------------

genes_a_evaluar <- unique(tabla_modelos$Gene)
total_genes <- length(genes_a_evaluar)
resultados_mmm <- list()

for (i in 1:total_genes) {
  g <- genes_a_evaluar[i]
  df_sub <- tabla_modelos %>% filter(Gene == g)
  n_muestras <- nrow(df_sub)
  
  if (i %% 50 == 0) {
    cat("Progreso:", i, "/", total_genes, "genes procesados (", round((i/total_genes)*100, 1), "%)\n")}
  
  if(nrow(df_sub) >= 30) { # Genes presentes en un mínimo de 30 muestras
    model_mix <- tryCatch({
      lmer(TPM ~ metilacion + as.factor(GrowthPattern) + MSIScore + PC1 + PC2 + (1 | OncotreeLineage), data = df_sub)
    }, error = function(e) { NULL })
    
    if(!is.null(model_mix)) {
      
      es_singular <- lme4::isSingular(model_mix) #Se evalúa la singularidad
      
      coefs <- broom.mixed::tidy(model_mix, effects = "fixed") %>%
        filter(term == "metilacion") #Se extrae el efecto de la variable predictora sobre la expresión
      
      aleatorios <- broom.mixed::tidy(model_mix, effects = "ran_pars") #Mide la variabilidad del linaje
      
      v_linaje   <- (aleatorios$estimate[aleatorios$group == "OncotreeLineage" & aleatorios$term == "sd__(Intercept)"])^2 #Varianza del linaje
      v_residual <- (aleatorios$estimate[aleatorios$group == "Residual" & aleatorios$term == "sd__Observation"])^2 #Varianza residual
      
      v_linaje   <- if (length(v_linaje) > 0 && !is.na(v_linaje)) v_linaje else 0 
      v_residual <- if (length(v_residual) > 0 && !is.na(v_residual)) v_residual else 0
      
      fijos_nombres <- names(lme4::fixef(model_mix))[-1] # Excluye el Intercepto
      v_fijos_lista <- numeric(length(fijos_nombres))
      names(v_fijos_lista) <- fijos_nombres
      
      for(fijo in fijos_nombres) {
        # Se calcula la varianza explicada por cada variable fija por separado
        v_fijos_lista[fijo] <- var(as.vector(lme4::getME(model_mix, "X")[, fijo] * lme4::fixef(model_mix)[fijo]))
      }
      
      v_meth   <- sum(v_fijos_lista[grepl("metilacion", names(v_fijos_lista))])
      v_growth <- sum(v_fijos_lista[grepl("GrowthPattern", names(v_fijos_lista))])
      v_msi    <- sum(v_fijos_lista[grepl("MSIScore", names(v_fijos_lista))])
      v_pc1    <- sum(v_fijos_lista[grepl("PC1", names(v_fijos_lista))])
      v_pc2    <- sum(v_fijos_lista[grepl("PC2", names(v_fijos_lista))])
      
      # Control de valores nulos o ausentes
      v_meth   <- if(is.na(v_meth)) 0 else v_meth
      v_growth <- if(is.na(v_growth)) 0 else v_growth
      v_msi    <- if(is.na(v_msi)) 0 else v_msi
      v_pc1    <- if(is.na(v_pc1)) 0 else v_pc1
      v_pc2    <- if(is.na(v_pc2)) 0 else v_pc2
      
      # Suma total (Efectos Fijos + Efecto Aleatorio + Varianza Residual)
      v_total_real <- v_meth + v_growth + v_msi + v_pc1 + v_pc2 + v_linaje + v_residual
      
      # Cálculo de la contribución porcentual
      if (v_total_real > 0) {
        pct_meth   <- (v_meth / v_total_real) * 100
        pct_growth <- (v_growth / v_total_real) * 100
        pct_msi    <- (v_msi / v_total_real) * 100
        pct_pc1    <- (v_pc1 / v_total_real) * 100
        pct_pc2    <- (v_pc2 / v_total_real) * 100
      } else { pct_meth <- 0; pct_growth <- 0; pct_msi <- 0; pct_pc1 <- 0; pct_pc2 <- 0 }
      
      # Se extraen los resultados
      if(nrow(coefs) > 0) {
        resultados_mmm[[g]] <- data.frame(
          Gene                         = g,
          estimate                     = coefs$estimate,
          statistic                    = coefs$statistic,
          p_value_metilacion           = coefs$p.value,
          Singular                     = es_singular,
          Pct_Var_Real_metilacion      = pct_meth,
          Pct_Var_Real_GrowthPattern   = pct_growth,
          Pct_Var_Real_MSIScore        = pct_msi,
          Pct_Var_Real_PC1             = pct_pc1,
          Pct_Var_Real_PC2             = pct_pc2,
          Pct_Var_Real_OncotreeLineage = (v_linaje / v_total_real) * 100,
          Pct_Var_Real_Residual        = (v_residual / v_total_real) * 100
        )}}}}

#Genes con muestras suficientes
genes_con_muestras_suficientes <- tabla_modelos %>% 
  group_by(Gene) %>% summarise(n = n(), .groups = "drop") %>% 
  filter(n >= 30) %>% # Genes con al menos 30 muestras
  nrow()

genes_modelados_exito <- length(resultados_mmm)# Genes resultantes del modelo

genes_con_error <- genes_con_muestras_suficientes - genes_modelados_exito # Genes que dieron error (tryCatch)

cat("Genes que tenían muestras suficientes de partida:", genes_con_muestras_suficientes, "\n") 
cat("Genes modelados con éxito total (sin errores):", genes_modelados_exito, "\n") 
cat("Total de genes que dieron error matemático en el lmer:", genes_con_error, "\n") 


# 7.2. Cálculo FDR 
#-------------------

df_resultados_unificados <- bind_rows(resultados_mmm) %>% filter(!is.na(Gene))

df_resultados_final <- df_resultados_unificados %>%
  mutate(ANOVA_P_metilacion_FDR = p.adjust(p_value_metilacion, method = "fdr")) %>%
  dplyr::select(Gene, estimate, statistic, ANOVA_P_metilacion_FDR, everything())


# 7.3. Genes que reducen su expresión mediante epigenética
#-----------------------------------------------------------

mmm_final <- df_resultados_final %>%
  filter(
    estimate < 0, # Correlación negativa
    ANOVA_P_metilacion_FDR < 0.05 # Solo genes significativos
  ) %>% arrange(estimate)

genes_silenciados_epigeneticamente <- unique(mmm_final$Gene)
cat("\nGenes que reducen su expresión mediante epigenética:", length(genes_silenciados_epigeneticamente), "\n") 

# Conteo de combinaciones Gen-Línea Celular con significancia estadística
combinaciones_reales_fdr <- tabla_modelos %>% filter(Gene %in% genes_silenciados_epigeneticamente) %>% nrow()

cat("Combinaciones tras FDR:", combinaciones_reales_fdr, "\n")  


# 7.4. Resumen de la varianza 
#-----------------------------

resumen_varianza <- mmm_final %>%
  dplyr::select(Gene, starts_with("Pct_Var_Real_")) %>%
  pivot_longer(cols = -Gene, names_to = "Variable", values_to = "Porcentaje") %>%
  group_by(Variable) %>%
  summarise(
    Varianza_Media_Explicada = mean(Porcentaje, na.rm = TRUE),
    Mediana_Varianza         = median(Porcentaje, na.rm = TRUE), 
    .groups = "drop"
  ) %>%
  mutate(
    Variable = case_when(
      Variable == "Pct_Var_Real_metilacion"     ~ "Metilación del Promotor (Epigenética)",
      Variable == "Pct_Var_Real_GrowthPattern"   ~ "Patrón de Crecimiento (GrowthPattern)",
      Variable == "Pct_Var_Real_MSIScore"        ~ "Inestabilidad de Microsatélites (MSI)",
      Variable == "Pct_Var_Real_PC1"             ~ "Componente Principal Epigenética 1 (PC1)",
      Variable == "Pct_Var_Real_PC2"             ~ "Componente Principal Epigenética 2 (PC2)",
      Variable == "Pct_Var_Real_OncotreeLineage" ~ "Linaje Tumoral / Tejido (Efecto Aleatorio / ICC)",
      Variable == "Pct_Var_Real_Residual"        ~ "Varianza Residual (No Explicada / Otros Factores)",
      TRUE                                       ~ as.character(Variable)
    )) %>% arrange(desc(Varianza_Media_Explicada))

print(as.data.frame(resumen_varianza))


#==================================================================================
# FASE 8: Genes con capacidad de silenciarse por hipermetilación de sus promotores 
#==================================================================================

genes_silenciados_epigeneticamente <- unique(mmm_final$Gene)

epi_sil_lines <- tabla_modelos %>%
  filter(Gene %in% genes_silenciados_epigeneticamente) %>%
  filter(metilacion >= 0.70 & TPM < umbral_matematico_exacto) %>% # Promotor hipermetilado y gen subexpresado 
  dplyr::select(ModelID, CCLEName, Gene, OncotreeLineage, GrowthPattern, MSIScore, PC1, PC2)

cat("Muestras moleculares con eventos de silenciamiento epigenético detectados:", nrow(epi_sil_lines), "\n") 

length(unique(epi_sil_lines$Gene)) #Genes silenciados epigenéticamente

# Determinación de la prevalencia por tejido
representacion_tejido <- epi_sil_lines %>% 
  group_by(Gene, OncotreeLineage) %>%
  summarise(Lineas_Silenciadas_En_Tejido = n(), .groups = "drop") %>%
  left_join(df_linajes_sincronizados %>% rename(Total_Lineas_Tejido = n), by = "OncotreeLineage") %>%
  mutate(Porcentaje_Tejido = (Lineas_Silenciadas_En_Tejido / Total_Lineas_Tejido) * 100)


# Número de tejidos en los que esta silenciado el gen
perfil_especificidad <- representacion_tejido %>%
  filter(Porcentaje_Tejido >= 10) %>% 
  group_by(Gene) %>%
  summarise(Total_Tejidos = n_distinct(OncotreeLineage), .groups = "drop")


#===========================================================
# FASE 9: Clasificación genes silenciados epigenéticamente 
#===========================================================

# Genes específicos de un tejido
genes_especificos <- perfil_especificidad %>%
  filter(Total_Tejidos == 1) %>% 
  pull(Gene)

n_tejidos_totales <- n_distinct(tabla_modelos$OncotreeLineage) # Número total de tipos de cáncer

cat("Número de genes específicos de tejido detectados:", length(genes_especificos), "\n")

# Genes comunes a más de la mitad de los tejidos
genes_universales <- perfil_especificidad %>%
  filter(Total_Tejidos >= (n_tejidos_totales / 2)) %>% 
  pull(Gene)

cat("Número de genes universales detectados:", length(genes_universales), "\n")


prevalencia <- epi_sil_lines %>%
  filter(Gene %in% unique(c(genes_especificos, genes_universales))) %>%
  group_by(Gene) %>%
  summarise(
    Muestras_Silenciadas = n_distinct(CCLEName), # Conteo de líneas celulares únicas
    Tejidos_Afectados    = n_distinct(OncotreeLineage), # Conteo de tejidos afectados
    .groups = "drop"
  ) %>% arrange(desc(Muestras_Silenciadas))

top_espe_205 <- prevalencia %>% filter(Gene %in% genes_especificos)

top_totales_538 <- prevalencia %>% filter(Gene %in% genes_universales)


# 9.1. Distribución de la carga epigenética 
# ------------------------------------------

genes_conj <- perfil_especificidad %>%  pull(Gene)

cell_line_especificos_profile <- epi_sil_lines %>%
  filter(Gene %in% genes_conj) %>%
  group_by(ModelID, OncotreeLineage) %>%
  summarise(n_genes_conj = n_distinct(Gene), .groups = "drop")

tabla_medianas_conj <- cell_line_especificos_profile %>%
  group_by(OncotreeLineage) %>%
  summarise(
    Mediana_Genes = median(n_genes_conj),
    Media_Genes   = round(mean(n_genes_conj), 2),
    Total_Muestras = n(),
    .groups = "drop"
  ) %>% arrange(desc(Mediana_Genes))

print(as.data.frame(tabla_medianas_conj))


plot_box_conj <- ggplot(cell_line_especificos_profile, aes(x = reorder(OncotreeLineage, n_genes_conj, FUN = median), y = n_genes_conj)) +
  geom_boxplot(fill = "#4A5568", alpha = 0.8, outlier.size = 1) +
  coord_flip() +
  scale_y_continuous(breaks = seq(0, 1200, by = 200)) + 
  theme_classic(base_size = 22) +
  labs(
    title = "Distribución de la carga epigenética específica por tejido",
    x = "Tejido de origen", y = "Genes específicos silenciados por muestra"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 24),
    axis.title = element_text(face = "bold", size = 22))

print(plot_box_conj) # Boxplot de genes silenciados epigenéticamente


# 9.2. Distribución de la carga epigenética de los genes específicos
# -------------------------------------------------------------------

cell_line_especificos_profile <- epi_sil_lines %>%
  filter(Gene %in% genes_especificos) %>%
  group_by(ModelID, OncotreeLineage) %>%
  summarise(n_genes_especificos = n_distinct(Gene), .groups = "drop")

tabla_medianas_espe <- cell_line_especificos_profile %>%
  group_by(OncotreeLineage) %>%
  summarise(
    Mediana_Genes = median(n_genes_especificos),
    Media_Genes   = round(mean(n_genes_especificos), 2),
    Total_Muestras = n(),
    .groups = "drop"
  ) %>% arrange(desc(Mediana_Genes))

print(as.data.frame(tabla_medianas_espe))

plot_box_espe <- ggplot(cell_line_especificos_profile, aes(x = reorder(OncotreeLineage, n_genes_especificos, FUN = median), y = n_genes_especificos)) +
  geom_boxplot(fill = "#4A5568", alpha = 0.8, outlier.size = 1) +
  coord_flip() + 
  theme_classic(base_size = 22) +
  labs(
    title = "Distribución de la carga epigenética específica por tejido (205 genes)",
    x = "Tejido de origen", y = "Genes específicos silenciados por muestra"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 24),
    axis.title = element_text(face = "bold", size = 22))

print(plot_box_espe) # Boxplot de Específicos


# 9.3. Distribución de la carga epigenética de los genes comunes
# ---------------------------------------------------------------

perfil_carga_538 <- epi_sil_lines %>% 
  filter(Gene %in% genes_universales) %>%
  group_by(ModelID, OncotreeLineage) %>%
  summarise(n_genes = n_distinct(Gene), .groups = "drop")

tabla_medianas_uni <- perfil_carga_538 %>%
  group_by(OncotreeLineage) %>%
  summarise(
    Mediana_Genes = median(n_genes),
    Media_Genes   = round(mean(n_genes), 2),
    Total_Muestras = n(),
    .groups = "drop"
  ) %>% arrange(desc(Mediana_Genes))

print(as.data.frame(tabla_medianas_uni))

plot_box_538 <- ggplot(perfil_carga_538, aes(x = reorder(OncotreeLineage, n_genes, FUN = median), y = n_genes)) +
  geom_boxplot(fill = "#4A5568", alpha = 0.8, outlier.size = 1) +
  coord_flip() + 
  theme_classic(base_size = 22) +
  labs(
    title = "Distribución de la carga universal inicial (538 genes)",
    x = "Tejido de origen", y = "Genes universales silenciados por muestra"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 24),
    axis.title = element_text(face = "bold", size = 22))

print(plot_box_538) # Boxplot de Universales Totales


#=======================================================
# FASE 10: Enriquecimiento funcional de genes comunes
#=======================================================

# Conversión de identificadores moleculares a ENTREZID

set.seed(42)

ids_univ_totales <- bitr(genes_universales, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)


# 10.1. Visualización de los procesos Gene Ontology en español
#---------------------------------------------------------------

# Análisis GO
set.seed(42)

go_univ_totales <- enrichGO(gene = ids_univ_totales$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)


if(!is.null(go_univ_totales) && nrow(as.data.frame(go_univ_totales)) > 0) {
  
  traduccion_universales <- c(
    "embryonic organ development"          = "Desarrollo de órganos embrionarios",
    "pattern specification process"        = "Proceso de especificación de patrones",
    "regionalization"                      = "Regionalización",
    "embryonic organ morphogenesis"        = "Morfogénesis de órganos embrionarios",
    "renal system development"             = "Desarrollo del sistema renal",
    "skeletal system morphogenesis"        = "Morfogénesis del sistema esquelético",
    "anterior/posterior pattern specification" = "Especificación del patrón anterior/posterior",
    "kidney development"                   = "Desarrollo del riñón",
    "embryonic skeletal system development" = "Desarrollo del sistema esquelético embrionario",
    "embryonic skeletal system morphogenesis" = "Morfogénesis del sistema esquelético embrionario")
  
  go_univ_esp <- go_univ_totales
  
  descripciones_originales_uni <- go_univ_esp@result$Description
  
  go_univ_esp@result$Description <- ifelse(
    descripciones_originales_uni %in% names(traduccion_universales),
    traduccion_universales[descripciones_originales_uni],
    descripciones_originales_uni)
  
  plot_go_universales <- dotplot(go_univ_esp, showCategory = 10) + 
    theme_classic(base_size = 18) +
    labs(
      title = "Procesos GO de los 538 genes comunes)",
      x = "Ratio génico",
      y = "Proceso biológico",
      color = "p.ajustado",
      size = "Conteo"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 24, hjust = 0.5),
      axis.title.x = element_text(face = "bold", size = 22, color = "black"),
      axis.title.y = element_text(face = "bold", size = 22, color = "black"),
      axis.text.y = element_text(size = 18, color = "black"),
      axis.text.x = element_text(size = 18, color = "black"),
      legend.title = element_text(face = "bold", size = 22) )
  
  print(plot_go_universales)
}


# 10.2. Análisis de tejidos y extracción de genes de los 10 GO principales
#--------------------------------------------------------------------------

if(!is.null(go_univ_esp) && nrow(as.data.frame(go_univ_esp)) > 0) {
  
  top_10_uni <- head(go_univ_esp@result$Description, 10) # Obtención de las 10 funciones principales
  
  df_base_genes <- go_univ_esp@result %>%
    filter(Description %in% top_10_uni) %>%
    dplyr::select(Description, geneID) %>%
    separate_rows(geneID, sep = "/") %>% 
    rename(Proceso_Biologico = Description, Gene = geneID) %>% 
    distinct()
  
  # Análisis de linajes y tejidos 
  tabla_modelos_clean <- tabla_modelos %>% 
    dplyr::select(Gene, OncotreeLineage) %>% 
    distinct()
  
  matriz_tejidos_go <- df_base_genes %>%
    inner_join(tabla_modelos_clean, by = "Gene", relationship = "many-to-many") %>% 
    group_by(Proceso_Biologico, OncotreeLineage) %>% 
    summarise(Num_Genes = n_distinct(Gene), .groups = "drop") %>%
    pivot_wider(names_from = OncotreeLineage, values_from = Num_Genes, values_fill = 0)

  print(as.data.frame(matriz_tejidos_go)) #conteo por tejido
  
  write.csv2(matriz_tejidos_go, "matriz_genes_por_tejido_go_limpio.csv", row.names = FALSE)
  
  # Extracción de los nombres de los genes
  tabla_lista_genes <- df_base_genes %>%
    group_by(Proceso_Biologico) %>%
    summarise(Genes_Asociados = paste(sort(unique(Gene)), collapse = "/"), .groups = "drop")
  
  print(as.data.frame(tabla_lista_genes)) #nombre de los genes por proceso
  
  write.csv2(tabla_lista_genes, "lista_genes_por_proceso_go.csv", row.names = FALSE)
  
} else {cat("No se puede realizar porque el objeto GO está vacío.\n")}


# 10.3. Visualización de rutas metabólicas (KEGG) en español
#------------------------------------------------------------

#Análisis de las rutas KEGG
set.seed(42)

kegg_univ_totales <- enrichKEGG(gene = ids_univ_totales$ENTREZID, organism = 'hsa', pAdjustMethod = "BH", pvalueCutoff = 0.05)


if(!is.null(kegg_univ_totales) && nrow(as.data.frame(kegg_univ_totales)) > 0) {
  
  traduccion_kegg_uni <- c(
    "Cytoskeleton in muscle cells"        = "Citoesqueleto en células musculares",
    "Integrin signaling"                   = "Señalización de integrinas",
    "Protein digestion and absorption"     = "Digestión y absorción de proteínas",
    "Leukocyte transendothelial migration" = "Migración transendotelial de leucocitos",
    "ECM-receptor interaction"             = "Interacción entre la ECM y el receptor",
    "Mineral absorption"                   = "Absorción de minerales")
  
  kegg_univ_esp <- kegg_univ_totales
  
  descripciones_originales_kegg <- kegg_univ_esp@result$Description
  
  kegg_univ_esp@result$Description <- ifelse(
    descripciones_originales_kegg %in% names(traduccion_kegg_uni),
    traduccion_kegg_uni[descripciones_originales_kegg],
    descripciones_originales_kegg)
  
  plot_kegg_universales <- dotplot(kegg_univ_esp, showCategory = 10) + 
    theme_classic(base_size = 18) +
    labs(
      title = "Rutas KEGG: Universales Totales Iniciales (538 Genes)",
      x = "Ratio génico",
      y = "Rutas KEGG",
      color = "p.ajustado",
      size = "Conteo"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 24, hjust = 0.5),
      axis.title.x = element_text(face = "bold", size = 22, color = "black"),
      axis.title.y = element_text(face = "bold", size = 22, color = "black"),
      axis.text.y = element_text(size = 18, color = "black"),
      axis.text.x = element_text(size = 18, color = "black"),
      legend.title = element_text(face = "bold", size = 22))
  
  print(plot_kegg_universales)
}


#=============================================================================================================
# FASE 11: Enriquecimiento funcional comparativo de los genes específicos linfoides según estadio madurativo
#=============================================================================================================

# Clasificación de líneas linfoides en estadio maduro e imaduro
mapeo_madurez_muestras <- model_global %>%
  dplyr::select(ModelID, OncotreePrimaryDisease) %>%
  distinct() %>%
  mutate(Estadio = case_when(
    # Estadio inmaduro
    OncotreePrimaryDisease %in% c("B-Cell Acute Lymphoblastic Leukemia", 
                                  "T-Lymphoblastic Leukemia/Lymphoma") ~ "Inmaduro",
    # Estadio maduro
    OncotreePrimaryDisease %in% c("Mature B-Cell Neoplasms", 
                                  "Mature T and NK Neoplasms", 
                                  "Non-Hodgkin Lymphoma", 
                                  "Hodgkin Lymphoma") ~ "Maduro",
    
    TRUE ~ NA_character_ # Exclusión de registros no clasificados
  )) %>% filter(!is.na(Estadio))


# Desglose líneas linfoides por estadio 
desglose_enfermedades <- mapeo_madurez_muestras %>%
  group_by(Estadio, OncotreePrimaryDisease) %>%
  summarise(
    Numero_Lineas_Celulares = n_distinct(ModelID),
    .groups = "drop"
  ) %>% arrange(Estadio, desc(Numero_Lineas_Celulares))

print(as.data.frame(desglose_enfermedades))


# Clasificación de los 188 genes específicos según su estadio
genes_clasificados <- epi_sil_lines %>%
  filter(Gene %in% genes_especificos) %>%
  inner_join(mapeo_madurez_muestras, by = "ModelID") %>%
  group_by(Gene, Estadio) %>%
  summarise(Muestras_Silenciadas = n(), .groups = "drop")

genes_linfoide_netos <- unique(genes_clasificados$Gene) # Vector 188 genes linfoides


# Se convierten los símbolos de los genes a ENTREZID
ids_espe_comparativo <- bitr(genes_linfoide_netos, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

genes_clasificados_entrez <- genes_clasificados %>%
  inner_join(ids_espe_comparativo, by = c("Gene" = "SYMBOL"))


# 11.1. Visualización comparativa de los procesos Gene Ontology en español
#--------------------------------------------------------------------------

# Análisis comparativo GO por estadio
set.seed(42)

go_comparativo <- compareCluster(
  ENTREZID ~ Estadio, 
  data = genes_clasificados_entrez, 
  fun = "enrichGO", 
  OrgDb = org.Hs.eg.db, 
  ont = "BP", 
  pAdjustMethod = "BH", 
  pvalueCutoff = 0.05, 
  readable = TRUE)

if(!is.null(go_comparativo) && nrow(as.data.frame(go_comparativo)) > 0) {
  
  go_comparativo_esp <- go_comparativo

  df_go_completo <- as.data.frame(go_comparativo_esp)
  
  top_10_reales <- df_go_completo %>% 
    arrange(p.adjust) %>% distinct(Description) %>%
    head(10) %>% pull(Description) # Se extraen los 10 procesos GO prncipales
  
  go_comparativo_esp@compareClusterResult <- go_comparativo_esp@compareClusterResult %>% 
    filter(Description %in% top_10_reales)
  
  traduccion_procesos <- c(
    "appendage morphogenesis"                         = "Morfogénesis de apéndices",
    "limb morphogenesis"                              = "Morfogénesis de extremidades",
    "embryonic skeletal system development"           = "Desarrollo del sistema esquelético embrionario",
    "cellular response to mineralocorticoid stimulus" = "Respuesta celular al estímulo mineralocorticoide",
    "appendage development"                           = "Desarrollo de apéndices",
    "limb development"                                = "Desarrollo de extremidades",
    "regulation of chondrocyte differentiation"       = "Regulación de la diferenciación de condrocitos",
    "morphogenesis of a branching structure"          = "Morfogénesis de una estructura ramificada",
    "cartilage development"                           = "Desarrollo del cartílago",
    "morphogenesis of a branching epithelium"         = "Morfogénesis de un epitelio ramificado")
  
  go_comparativo_esp@compareClusterResult$Description <- traduccion_procesos[go_comparativo_esp@compareClusterResult$Description]
  
  plot_comparativo <- dotplot(
    go_comparativo_esp, 
    showCategory = 10, # Límite de 10 categorías
    by = "count"
  ) + 
    theme_classic(base_size = 18) +
    labs(
      title = "Procesos GO comparativos de los tumores linfoides",
      x = "Estadio madurativo celular",
      y = "Proceso biológico",
      color = "p.ajustado",
      size = "Conteo"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 24, hjust = 0.5),
      axis.title.x = element_text(face = "bold", size = 22, color = "black"),
      axis.title.y = element_text(face = "bold", size = 22, color = "black"),
      axis.text.y = element_text(face = "plain", size = 18, color = "black"),
      axis.text.x = element_text(face = "plain", size = 18, color = "black"),
      legend.title = element_text(face = "bold", size = 22))
  
  print(plot_comparativo)
  
  # Extracción de genes por proceso GO y estadio madurativo
  df_go_filtrado <- df_go_completo %>%
    filter(Description %in% top_10_reales) %>%
    dplyr::select(Cluster, ID, Description, geneID) %>%
    rename(Estadio = Cluster) %>% 
    separate_rows(geneID, sep = "/") %>% # Se separan los genes por filas 
    rename(Gene = geneID) %>%
    distinct()
  
  genes_por_proceso_y_estadio <- df_go_filtrado %>%
    group_by(Description, Estadio) %>% # Se agrupa por proceso biológico y estadio 
    summarise(
      Total_Genes = n(),
      Lista_Genes = paste(sort(Gene), collapse = ", "),
      .groups = "drop"
    ) %>%
    arrange(Description, Estadio)
  
  print(as.data.frame(genes_por_proceso_y_estadio))
  
  write.csv2(genes_por_proceso_y_estadio, "top10_go_genes_por_estadio.csv", row.names = FALSE)
  
} else {cat("No se encontraron términos significativos en el análisis de compareCluster.\n")}


# 11.2. Visualización comparativa de los procesos Gene Ontology en español
#--------------------------------------------------------------------------

# Análisis comparativo de las rutas KEGG
set.seed(42)

kegg_comparativo <- compareCluster(
  ENTREZID ~ Estadio, 
  data = genes_clasificados_entrez, 
  fun = "enrichKEGG", 
  organism = "hsa", # Código para humanos en KEGG
  pAdjustMethod = "BH", 
  pvalueCutoff = 0.05)


if(!is.null(kegg_comparativo) && nrow(as.data.frame(kegg_comparativo)) > 0) {
  
  traduccion_kegg <- c(
    "GnRH secretion" = "Secreción de GnRH",
    "Parathyroid hormone synthesis, secretion and action" = "Síntesis, secreción y acción de la hormona paratiroidea")
  
  kegg_comparativo_esp <- kegg_comparativo
  
  kegg_comparativo_esp@compareClusterResult$Description <- traduccion_kegg[kegg_comparativo_esp@compareClusterResult$Description]
  
  kegg_comparativo_esp@compareClusterResult$p.adjust <- round(kegg_comparativo_esp@compareClusterResult$p.adjust, 3)
  
  plot_kegg_comp <- dotplot(kegg_comparativo_esp, showCategory = 10) + 
    theme_classic(base_size = 18) + 
    scale_x_discrete(labels = c("Maduro")) +
    labs(
      title = "Rutas KEGG comparativas del linaje linfoide",
      x = "Estadio madurativo celular", 
      y = "Rutas KEGG",
      color = "p.ajustado",
      size = "Ratio génico"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 24, hjust = 0.5),
      axis.title.x = element_text(face = "bold", size = 22, color = "black"),
      axis.title.y = element_text(face = "bold", size = 22, color = "black"),
      axis.text.y = element_text(face = "plain", size = 22, color = "black"),
      axis.text.x = element_text(face = "plain", size = 22, color = "black"),
      legend.title = element_text(face = "bold", size = 22),
      legend.text = element_text(size = 18, color = "black"))
  
  print(plot_kegg_comp)
  
} else {cat("No se encontraron rutas KEGG enriquecidas significativamente con este p-value.\n")}


#==============================================
# FASE 12: Viabilidad celular ante un fármaco 
#==============================================

# 12.1. Cribado de compuestos por representatividad y variabilidad 
#------------------------------------------------------------------

columnas_farma_num <- intersect(colnames(farma1), unique(tabla_modelos$ModelID))

farmacos_variables <- farma1 %>%
  # El fármaco debe poseer registros de viabilidad en al menos el 50% de la cohorte
  filter(rowSums(is.na(.[, columnas_farma_num])) < length(columnas_farma_num) * 0.5) %>% 
  # Cálculo de la desviación estándar de la respuesta celular por cada compuesto
  mutate(SD_Viability = apply(.[, columnas_farma_num, drop = FALSE], 1, sd, na.rm = TRUE)) %>%
  filter(SD_Viability >= 0.15) %>% # Filtro de varianza
  pull(IDs)

farma1_filtrado_viabilidad <- farma1 %>% filter(IDs %in% farmacos_variables)


# 12.2. Filtrado para obtener solo fármacos que tengan MOA o Target
#--------------------------------------------------------------------

farma2_filtrado <- farma2 %>%
  filter((!is.na(repurposing_target) & repurposing_target != "") | (!is.na(MOA) & MOA != ""))


farma_final <- farma1_filtrado_viabilidad %>% 
  semi_join(farma2_filtrado, by = "IDs") %>%
  dplyr::select(IDs, all_of(columnas_farma_num)) %>%
  pivot_longer(
    cols = -IDs,
    names_to = "ModelID",
    values_to = "viability"
  ) %>%
  filter(!is.na(viability)) %>%
  left_join(distinct(farma2_filtrado, IDs, .keep_all = TRUE), by = "IDs")


#12.3. Reporte
#--------------

cat("Fármacos iniciales:", length(unique(farma1$IDs)), "\n")
cat("Fármacos tras cribado:", length(unique(farma1_filtrado_viabilidad$IDs)), "\n") 
cat("Fármacos finales tras el cribado:", length(unique(farma_final$IDs)), "\n")


#===================================================
# FASE 13: Modelo lineal multivariante regularizado
#===================================================

# 13.1. Configuración del cluster multinúcleo
#----------------------------------------------

try(parallel::stopCluster(cl), silent = TRUE) #Detiene el cluster previo si estaba abierto

closeAllConnections() #Cierra las conexiones abiertas

gc(verbose = FALSE) #Libera memoria RAM

num_nucleos <- parallel::detectCores() - 2 # 2 núcleos libres

cl <- parallel::makeCluster(num_nucleos) #Se crea el cluster

doParallel::registerDoParallel(cl) #Bucles paralelos

doRNG::registerDoRNG(42) #Se fija semilla


# 13.2. Optimización pre-bucle
# -----------------------------

todos_los_farmacos <- unique(farma_final$IDs)

total_f <- length(todos_los_farmacos)


# Estructura clínica basal 
estructura_clinica_base <- tabla_modelos %>% 
  distinct(ModelID, .keep_all = TRUE) %>% 
  dplyr::select(ModelID, OncotreeLineage, GrowthPattern, MSIScore, PC1, PC2) %>%
  left_join(mapeo_madurez_muestras %>% dplyr::select(ModelID, Estadio), by = "ModelID") %>%
  mutate(
    OncotreeLineage = case_when(
      !is.na(Estadio) & Estadio == "Inmaduro" ~ "Linfoide_Inmaduro",
      !is.na(Estadio) & Estadio == "Maduro"   ~ "Linfoide_Maduro",
      TRUE ~ as.character(OncotreeLineage)
    ),
    OncotreeLineage = as.factor(OncotreeLineage),
    GrowthPattern = if_else(as.character(GrowthPattern) %in% c("1", "Adherent", "Mixed"), 1, 0) )


# Se identifican qué muestras pertenecen a cada estadio linfoide
muestras_inmaduras <- mapeo_madurez_muestras %>% filter(Estadio == "Inmaduro") %>% pull(ModelID)
muestras_maduras   <- mapeo_madurez_muestras %>% filter(Estadio == "Maduro") %>% pull(ModelID)


# Matriz ómica basal
matriz_predictores_lasso_base <- methyl_final %>%
  filter(Gene %in% perfil_especificidad$Gene) %>% 
  mutate(metilacion = as.numeric(as.character(metilacion))) %>%
  filter(!is.na(metilacion)) %>%
  pivot_wider(id_cols = ModelID, names_from = Gene, values_from = metilacion, values_fn = mean)


# Matrices independientes para cada estadio linfoide
matriz_interaccion_inmaduro <- matriz_predictores_lasso_base
matriz_interaccion_maduro   <- matriz_predictores_lasso_base

# Se ponen a 0 la metilación si la muestra no es inmadura
matriz_interaccion_inmaduro[!(matriz_interaccion_inmaduro$ModelID %in% muestras_inmaduras), -1] <- 0
colnames(matriz_interaccion_inmaduro)[-1] <- paste0(colnames(matriz_interaccion_inmaduro)[-1], "_x_LinfoideInmaduro") #Renombra los genes

# Se ponen a 0 la metilación si la muestra no es maduras
matriz_interaccion_maduro[!(matriz_interaccion_maduro$ModelID %in% muestras_maduras), -1] <- 0
colnames(matriz_interaccion_maduro)[-1] <- paste0(colnames(matriz_interaccion_maduro)[-1], "_x_LinfoideMaduro") #Renombra los genes


# Matriz total de predictores
matriz_predictores_lasso <- matriz_predictores_lasso_base %>%
  inner_join(matriz_interaccion_inmaduro, by = "ModelID") %>%
  inner_join(matriz_interaccion_maduro, by = "ModelID")


# Lista de genes candidatos
genes_candidatos_todos <- intersect(unique(perfil_especificidad$Gene), colnames(matriz_predictores_lasso_base))


# Matriz de diseño clínica dummificada 
X_clinicas_matriz_global <- model.matrix(~ GrowthPattern + MSIScore + PC1 + PC2 + OncotreeLineage, 
                                         data = estructura_clinica_base)

rownames(X_clinicas_matriz_global) <- estructura_clinica_base$ModelID

X_clinicas_matriz_global <- X_clinicas_matriz_global[, -1] # Se quita el intercepto


# Indexación de la tabla de fármacos 
farma_final_preparado <- farma_final %>%
  mutate(viability = as.numeric(as.character(viability))) %>%
  filter(!is.na(viability)) %>%
  dplyr::select(IDs, ModelID, viability)

lista_farmacos_split <- split(farma_final_preparado, farma_final_preparado$IDs)


# Variables para la paralelización
mis_paquetes <- c("glmnet", "dplyr", "tidyr")

mis_variables <- c("todos_los_farmacos", "lista_farmacos_split", "estructura_clinica_base", 
                   "matriz_predictores_lasso", "genes_candidatos_todos", "X_clinicas_matriz_global",
                   "genes_linfoide_netos")


# 13.3. Bucle paralelizado 
# -------------------------

tabla_resistencia_global <- foreach(
  i         = 1:total_f, 
  .combine  = bind_rows,
  .packages = mis_paquetes,
  .export   = mis_variables,
  .options.RNG = 42  
) %dopar% {
  
  #Se identifica el fármaco actual y se extraen sus datos
  farmaco_objetivo <- todos_los_farmacos[i]
  datos_farmaco_sub <- lista_farmacos_split[[farmaco_objetivo]]
  
  #Si el fármaco no tiene datos disponibles, se pasa al siguiente
  if (is.null(datos_farmaco_sub) || nrow(datos_farmaco_sub) == 0) return(NULL)
  
  # Se alinean los datos del fármaco actual
  tabla_lasso_final <- estructura_clinica_base %>% 
    inner_join(datos_farmaco_sub, by = "ModelID") %>%
    inner_join(matriz_predictores_lasso, by = "ModelID")
  
  # Se aíslan los datos clínicos del fármaco actual
  muestras_efectivas <- tabla_lasso_final$ModelID
  Y_vector <- as.matrix(tabla_lasso_final$viability)
  X_base <- X_clinicas_matriz_global[muestras_efectivas, , drop = FALSE]
  
  # Se eliminan las variables clínicas sin varianza para el fármaco actual
  columnas_validas <- apply(X_base, 2, function(x) length(unique(x)) > 1)
  X_base <- X_base[, columnas_validas, drop = FALSE]
  
  resultados_screening_mixto <- numeric(length(genes_candidatos_todos))
  names(resultados_screening_mixto) <- genes_candidatos_todos
  
  # Regresión lineal gen a gen 
  for (g in genes_candidatos_todos) { 
    gene_vector <- as.matrix(tabla_lasso_final[, g, drop = FALSE])
    
    #Imputa valores faltantes con la mediana del gen
    if (any(is.na(gene_vector))) {
      gene_vector[is.na(gene_vector)] <- median(gene_vector, na.rm = TRUE)}
    
    #Se descarta el gen si o tiene varianza
    if (var(gene_vector, na.rm = TRUE) == 0) { 
      resultados_screening_mixto[g] <- 0 
      next} 
    
    #Ajusta regresión lineal por mínimos cuadrados ordinarios
    X_completa_lm <- cbind(gene_vector, X_base)
    fit_g <- lm.fit(x = X_completa_lm, y = Y_vector)
    
    p <- fit_g$rank 
    
    if (p > 0) { 
      # Calcula el error estándar del modelo
      se_g <- tryCatch({
        resid_var <- sum(fit_g$residuals^2) / (length(Y_vector) - p)
        X_inv_diag <- diag(solve(t(X_completa_lm) %*% X_completa_lm))
        sqrt(resid_var * X_inv_diag)
      }, error = function(e) NULL)
      
      if (!is.null(se_g)) {
        # Se extrae el estadístico t asociado al gen en valor absoluto
        resultados_screening_mixto[g] <- abs(fit_g$coefficients / se_g)[1] 
      } else { resultados_screening_mixto[g] <- 0 } 
    } else { resultados_screening_mixto[g] <- 0 }
  }
  
  # Se selecciona el top 15 de genes con mayor asociación (estadístico t más alto)
  genes_validados_etapa1 <- names(sort(resultados_screening_mixto, decreasing = TRUE)[1:min(15, length(resultados_screening_mixto))])
  
  # Si no hay genes válidos, finaliza la iteración
  if(length(genes_validados_etapa1) == 0) return(NULL)
  
  # Se filtran las filas clínicas correspondientes a las muestras efectivas
  X_clinicas_matriz <- X_clinicas_matriz_global[muestras_efectivas, , drop = FALSE]
  
  # Se eliminan variables clínicas que no tengan varianza
  X_clinicas_matriz <- X_clinicas_matriz[, apply(X_clinicas_matriz, 2, function(x) length(unique(x)) > 1), drop = FALSE]
  
  variables_lasso_bloque <- c(
    genes_validados_etapa1,
    paste0(intersect(genes_validados_etapa1, genes_linfoide_netos), "_x_LinfoideInmaduro"),
    paste0(intersect(genes_validados_etapa1, genes_linfoide_netos), "_x_LinfoideMaduro")
  )
  variables_lasso_bloque <- intersect(variables_lasso_bloque, colnames(tabla_lasso_final))
  
  # Se extrae la submatriz ómica y de interacciones para el modelo
  matriz_genes_limpia <- tabla_lasso_final[, variables_lasso_bloque, drop = FALSE]
  
  # Imputación de NAs con la mediana en la matriz ómica y de interacciones
  for (col in colnames(matriz_genes_limpia)) {
    if (any(is.na(matriz_genes_limpia[, col]))) { 
      mediana_val <- median(matriz_genes_limpia[, col], na.rm = TRUE)
      matriz_genes_limpia[, col] <- ifelse(is.na(matriz_genes_limpia[, col]), mediana_val, matriz_genes_limpia[, col])
    }
  }
  X_genes_limpios <- as.matrix(matriz_genes_limpia)
  X_mat_completa <- cbind(X_clinicas_matriz, X_genes_limpios)
 
  # Se penaliza las ómicas con 1 y las clínicas con 0
  vector_penalizacion <- c(rep(0, ncol(X_clinicas_matriz)), rep(1, ncol(X_genes_limpios)))
  
  set.seed(123) 
  n_muestras_efectivas <- nrow(X_mat_completa)
  # Se asignan las muestras a 5 grupos balanceados
  mis_folds_fijos <- sample(rep(1:5, length.out = n_muestras_efectivas))
  
  # Se ejecuta la regresión Elastic Net
  cv_lasso <- tryCatch(
    cv.glmnet(x = X_mat_completa, y = Y_vector, alpha = 0.5, 
              penalty.factor = vector_penalizacion, standardize = TRUE, 
              foldid = mis_folds_fijos), 
    error = function(e) NULL)
  
  if (!is.null(cv_lasso)) { # Si el modelo converge, se extraen y guardan los resultados
    # Se obtienen los coeficientes asociados al lambda óptimo
    coeficientes_mat <- as.matrix(coef(cv_lasso, s = "lambda.min"))
    coeficientes_genes <- coeficientes_mat[intersect(rownames(coeficientes_mat), 
                                                     variables_lasso_bloque), , drop = FALSE]
    
    data.frame(
      Farmaco     = farmaco_objetivo,
      Gen         = rownames(coeficientes_genes),
      Coeficiente = as.numeric(coeficientes_genes[, 1]))
  } else { return(NULL) }
}


# 13.4. Apagado del cluster
# --------------------------

parallel::stopCluster(cl) # Se apaga el clúster en paralelo y se liberan los núcleos del procesador

rm(lista_farmacos_split, X_clinicas_matriz_global) #Se eliminan objetos para liberar memoria

gc(verbose = FALSE) #Libera memoria RAM


# ===================================================
# FASE 14: Extracción de datos del linaje linfoide
# ===================================================

# Homogeneización del diccionario
diccionario_completo <- farma2_filtrado %>%
  distinct(IDs, .keep_all = TRUE) %>%
  mutate(IDs_Limpio = gsub("^BRD:", "", as.character(IDs))) %>% 
  dplyr::select(IDs_Limpio, Nombre_Real_Farmaco = Drug.Name, Mecanismo_Clinico = MOA, Diana_Genetica = repurposing_target)


# Unión con metadatos 
tabla_resistencia_mapeada <- tabla_resistencia_global %>%
  mutate(Farmaco_Limpio = gsub("^BRD:", "", as.character(Farmaco))) %>% 
  left_join(diccionario_completo, by = c("Farmaco_Limpio" = "IDs_Limpio")) %>%
  mutate(
    Farmaco_Final = if_else(is.na(Nombre_Real_Farmaco) | Nombre_Real_Farmaco == "", as.character(Farmaco_Limpio), as.character(Nombre_Real_Farmaco)),
    MOA_Final     = if_else(is.na(Mecanismo_Clinico) | Mecanismo_Clinico == "", "Experimental / Desconocido", as.character(Mecanismo_Clinico)),
    Target_Final  = if_else(is.na(Diana_Genetica) | Diana_Genetica == "", "No anotado", as.character(Diana_Genetica))
  ) %>%
  dplyr::select(Gen, Farmaco = Farmaco_Final, Coeficiente, MOA = MOA_Final, Target = Target_Final)


# Se extraen los coeficientes del modelo regularizado
coeficientes_activos_modelo <- tabla_resistencia_mapeada %>%
  filter(!is.na(Coeficiente) & Coeficiente != 0) %>% 
  mutate(
    Estadio = case_when(
      grepl("_x_LinfoideInmaduro", Gen) ~ "Inmaduro",
      grepl("_x_LinfoideMaduro", Gen)   ~ "Maduro",
      TRUE ~ "Otro"
    ),
    Nombre_Gen = gsub("_x_LinfoideInmaduro|_x_LinfoideMaduro", "", Gen)
  ) %>%
  filter(Estadio != "Otro") %>%
  dplyr::select(Estadio, Gen = Nombre_Gen, Farmaco, Coeficiente, MOA) %>%
  distinct()


perfil_muestras_silenciadas_reales <- tabla_modelos %>%
  filter(Gene %in% coeficientes_activos_modelo$Gen) %>%
  dplyr::select(ModelID, Gen = Gene, metilacion, TPM) %>%
  mutate(
    metilacion = as.numeric(as.character(metilacion)),
    TPM = as.numeric(as.character(TPM))
  ) %>%
  filter(!is.na(metilacion) & !is.na(TPM)) %>%
  inner_join(mapeo_madurez_muestras %>% dplyr::select(ModelID, Estadio_Muestra = Estadio), by = "ModelID") %>%
  filter(metilacion >= 0.70 & TPM < umbral_matematico_exacto) %>% #Hipermetilados y silenciados
  group_by(Estadio = Estadio_Muestra, Gen) %>%
  summarise(Muestras_Con_Silenciamiento_Puro = n(), .groups = "drop")


farma_linfoides_filtrado <- coeficientes_activos_modelo %>%
  inner_join(perfil_muestras_silenciadas_reales, by = c("Estadio", "Gen")) %>%
  filter(Muestras_Con_Silenciamiento_Puro > 0)


# Exclusión de mecanismos fisiológicos o no-oncológicos generales
moa_excluir <- c("ADRENERGIC", "SEROTONIN", "DIURETIC", "NEUROMUSCULAR", "BLOOD PRESSURE", 
                 "PSYCHOTIC", "ANESTHETIC", "OPIOID", "CELL WALL", "BACTERIA", "BACTERIAL", 
                 "ANTIFUNGAL", "ANTIOXIDANT", "ANTIVIRAL", "CALCIUM CHANNEL", "GLUTAMATE", 
                 "HISTAMINE", "DOPAMINE", "ACETYLCHOLINE", "SECRETAC", "GROWTH HORMONE", 
                 "NEPRILYSIN", "CANNABINOID", "SENSING", "UROTENSIN", "PROGESTERONE", 
                 "CCK", "ANDROGEN")

farma_linfoides_filtrado <- farma_linfoides_filtrado %>%
  filter(!str_detect(MOA, regex(str_c(moa_excluir, collapse = "|"), ignore_case = TRUE))) %>%
  mutate(Estadio = if_else(Estadio == "Inmaduro", "Estadio linfoide inmaduro", "Estadio linfoide maduro"))


# Se aíslan los 5 MOA más fuertes por estadio
mecanismos_ganadores_reales <- farma_linfoides_filtrado %>% 
  group_by(Estadio, MOA) %>% 
  summarise(Fuerza_Efecto_Media = mean(abs(Coeficiente), na.rm = TRUE), .groups = "drop_last") %>%
  arrange(desc(Fuerza_Efecto_Media), MOA, .by_group = TRUE) %>% 
  slice_head(n = 5) %>% 
  ungroup() %>% 
  dplyr::select(Estadio, MOA)


tabla_frecuencias_reales <- farma_linfoides_filtrado %>% 
  semi_join(mecanismos_ganadores_reales, by = c("Estadio", "MOA")) %>%
  mutate(Efecto = if_else(Coeficiente < 0, "Metilación induce sensibilidad", "Metilación induce resistencia")) %>%
  group_by(Estadio, MOA, Efecto) %>% 
  summarise(Genes_Unicos = n_distinct(Gen), Farmacos_Unicos = n_distinct(Farmaco), 
            Metrica_Normalizada = Genes_Unicos / Farmacos_Unicos, .groups = "drop")


# Traducción 
objeto_grafico_linfoides <- tabla_frecuencias_reales %>%
  mutate(MOA_Esp = case_when(
    str_detect(MOA, regex("MCL1 INHIBITOR", ignore_case = TRUE))          ~ "Inhibidor de MCL1",
    str_detect(MOA, regex("EXPORTIN", ignore_case = TRUE))                ~ "Inhibidor de exportina",
    str_detect(MOA, regex("PET RADIOTRACER", ignore_case = TRUE))         ~ "Radiotrazador para PET",
    str_detect(MOA, regex("SODIUM/HYDROGEN", ignore_case = TRUE))         ~ "Inhibidor de Na+/H+ y canales TRP",
    str_detect(MOA, regex("PPAR RECEPTOR INVERSE", ignore_case = TRUE))   ~ "Agonista inverso del receptor PPAR",
    str_detect(MOA, regex("MITOCHONDRIAL COMPLEX I", ignore_case = TRUE)) ~ "Inhibidor del complejo mitocondrial I",
    str_detect(MOA, regex("MAP KINASE ACTIVATOR", ignore_case = TRUE))    ~ "Activador de la vía MAP quinasa",
    str_detect(MOA, regex("SREBP", ignore_case = TRUE))                   ~ "Inhibidor de SREBP",
    str_detect(MOA, regex("AURORA KINASE", ignore_case = TRUE))           ~ "Inhibidor de Aurora quinasa",
    str_detect(MOA, regex("IKK INHIBITOR", ignore_case = TRUE))           ~ "Inhibidor de IKK",
    str_detect(MOA, regex("RHO ASSOCIATED KINASE", ignore_case = TRUE))   ~ "Inhibidor de quinasa asociada a Rho",
    str_detect(MOA, regex("MICROTUBULE INHIBITOR", ignore_case = TRUE))   ~ "Inhibidor de microtúbulos",
    str_detect(MOA, regex("MDM INHIBITOR", ignore_case = TRUE))           ~ "Inhibidor de MDM",
    str_detect(MOA, regex("SIGMA RECEPTOR LIGAND", ignore_case = TRUE))   ~ "Ligando del receptor Sigma",
    str_detect(MOA, regex("JAK INHIBITOR", ignore_case = TRUE))           ~ "Inhibidor de JAK",
    str_detect(MOA, regex("RAF INHIBITOR", ignore_case = TRUE))           ~ "Inhibidor de RAF",
    str_detect(MOA, regex("MYELOPEROXIDASE", ignore_case = TRUE))         ~ "Inhibidor de mieloperoxidasa",
    str_detect(MOA, regex("MONOCARBOXYLATE", ignore_case = TRUE))         ~ "Inhibidor de transportadores de monocarboxilato",
    TRUE ~ as.character(MOA) 
  )) %>%
  mutate(Frecuencia_Grafico = if_else(Efecto == "Metilación induce resistencia", Metrica_Normalizada, -Metrica_Normalizada)) %>%
  group_by(Estadio, MOA_Esp) %>% 
  mutate(Total_Orden = sum(abs(Frecuencia_Grafico))) %>% 
  ungroup()


#Visualización
grafico_linfoides_comparativo <- ggplot(objeto_grafico_linfoides, aes(x = reorder(MOA_Esp, Total_Orden), y = Frecuencia_Grafico, fill = Efecto)) +
  geom_bar(stat = "identity", alpha = 0.85, width = 0.65) + 
  coord_flip() + 
  facet_wrap(~Estadio, scales = "free", ncol = 2) +
  scale_fill_manual(values = c("Metilación induce sensibilidad" = "#2B6CB0", "Metilación induce resistencia" = "#C53030")) +
  scale_y_continuous(labels = abs, limits = c(-1.5, 1.5)) + 
  theme_minimal(base_size = 14) + 
  labs(title = "Mecanismos de acción farmacológica determinantes por estadio linfoide", 
       x = "Mecanismo de acción clínico (MOA)", y = "Densidad de interacciones (Genes predictores por fármaco único)", fill = NULL) +
  theme(plot.title = element_text(face = "bold", size = 19, color = "#1A202C"), axis.title = element_text(face = "bold", size = 19), 
        axis.text.y = element_text(size = 17, color = "#2D3748"), axis.text.x = element_text(size = 17, color = "#2D3748"), 
        strip.text = element_text(face = "bold", size = 17, color = "#1A202C"), legend.text = element_text(size = 17), 
        legend.position = "bottom", panel.grid.minor = element_blank())

print(grafico_linfoides_comparativo)


# 14.1. REPORTE 
# ---------------

# Traducción
farma_linfoides_traducido <- farma_linfoides_filtrado %>%
  mutate(Mecanismo_Grafico = case_when(
    str_detect(MOA, regex("MCL1 INHIBITOR", ignore_case = TRUE))          ~ "Inhibidor de MCL1",
    str_detect(MOA, regex("EXPORTIN", ignore_case = TRUE))                ~ "Inhibidor de exportina",
    str_detect(MOA, regex("PET RADIOTRACER", ignore_case = TRUE))         ~ "Radiotrazador para PET",
    str_detect(MOA, regex("SODIUM/HYDROGEN", ignore_case = TRUE))         ~ "Inhibidor de Na+/H+ y canales TRP",
    str_detect(MOA, regex("MAP KINASE ACTIVATOR", ignore_case = TRUE))    ~ "Activador de la vía MAP quinasa",
    str_detect(MOA, regex("SREBP", ignore_case = TRUE))                   ~ "Inhibidor de SREBP",
    str_detect(MOA, regex("AURORA KINASE", ignore_case = TRUE))           ~ "Inhibidor de Aurora quinasa",
    str_detect(MOA, regex("IKK INHIBITOR", ignore_case = TRUE))           ~ "Inhibidor de IKK",
    str_detect(MOA, regex("RHO ASSOCIATED KINASE", ignore_case = TRUE))   ~ "Inhibidor de quinasa asociada a Rho",
    str_detect(MOA, regex("MICROTUBULE INHIBITOR", ignore_case = TRUE))   ~ "Inhibidor de microtúbulos",
    TRUE ~ as.character(MOA) 
  ))


referencia_grafica_limpia <- objeto_grafico_linfoides %>%
  dplyr::select(Estadio, Mecanismo_Grafico = MOA_Esp, Efecto_Grafico = Efecto) %>%
  distinct()


# Extracción de metadatos ómicos para las muestras que pasan el filtro
metadatos_omicos_filtrados_reales <- tabla_modelos %>%
  filter(Gene %in% farma_linfoides_traducido$Gen) %>%
  dplyr::select(ModelID, Gen = Gene, metilacion, TPM) %>%
  mutate(
    metilacion = as.numeric(as.character(metilacion)),
    TPM = as.numeric(as.character(TPM))
  ) %>%
  filter(!is.na(metilacion) & !is.na(TPM)) %>%
  inner_join(mapeo_madurez_muestras %>% dplyr::select(ModelID, Estadio_Muestra = Estadio), by = "ModelID") %>%
  mutate(Estadio = if_else(Estadio_Muestra == "Inmaduro", "Estadio linfoide inmaduro", "Estadio linfoide maduro")) %>%
  filter(metilacion >= 0.70 & TPM < umbral_matematico_exacto) %>% #Hipermetilación y silenciamiento
  group_by(Estadio, Gen) %>%
  summarise(
    Media_Metilacion = mean(metilacion, na.rm = TRUE),
    Media_Expresion_TPM = mean(TPM, na.rm = TRUE),
    Mediana_Metilacion = median(metilacion, na.rm = TRUE),
    Mediana_Expresion_TPM = median(TPM, na.rm = TRUE),
    Lineas_Celulares_Silenciadas = n(), 
    .groups = "drop"
  )


reporte_tfm_linfoides_con_coef <- farma_linfoides_traducido %>%
  dplyr::select(Estadio, Nombre_Gen = Gen, Farmaco, Coeficiente, Mecanismo_Grafico) %>%
  # Se asigna el efecto clínico según el signo del modelo
  mutate(Efecto_Real_Signo = if_else(Coeficiente < 0, "Metilación induce sensibilidad", "Metilación induce resistencia")) %>%
  inner_join(
    referencia_grafica_limpia, 
    by = c("Estadio", "Mecanismo_Grafico", "Efecto_Real_Signo" = "Efecto_Grafico")
  ) %>%
  left_join(metadatos_omicos_filtrados_reales, by = c("Estadio", "Nombre_Gen" = "Gen")) %>%
  mutate(Efecto_Clinico = Efecto_Real_Signo) %>%
  dplyr::select(-Efecto_Real_Signo) %>%
  distinct(Estadio, Mecanismo_Grafico, Nombre_Gen, Farmaco, Efecto_Clinico, .keep_all = TRUE) %>%
  group_by(Estadio, Mecanismo_Grafico, Efecto_Clinico) %>%
  arrange(Estadio, Mecanismo_Grafico, Efecto_Clinico, desc(abs(Coeficiente))) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  arrange(desc(Estadio), desc(abs(Coeficiente))) %>%
  dplyr::select(
    Estadio, Mecanismo_Grafico, Nombre_Gen, Farmaco, Coeficiente, Efecto_Clinico,
    Media_Metilacion, Media_Expresion_TPM, Mediana_Metilacion, Mediana_Expresion_TPM,
    Lineas_Celulares_Silenciadas
  )


write.csv2(reporte_tfm_linfoides_con_coef, "Reporte_TFM_Linfoides_farma.csv", row.names = FALSE)

print(as.data.frame(reporte_tfm_linfoides_con_coef))


# ====================================================================
# FASE 15: Extracción de datos comunes a todos los tumores analizados
# ====================================================================

umbral_universal  <- n_tejidos_totales / 2
genes_universales_perfil <- perfil_especificidad %>% filter(Total_Tejidos >= umbral_universal) %>% pull(Gene)

# Se extraen los coeficientes del modelo regularizado
coeficientes_activos_universales <- tabla_resistencia_mapeada %>%
  filter(!is.na(Coeficiente) & Coeficiente != 0) %>% 
  filter(Gen %in% genes_universales_perfil) %>%
  dplyr::select(Gen, Farmaco, Coeficiente, MOA) %>%
  distinct()


perfil_muestras_universales_reales <- tabla_modelos %>%
  filter(Gene %in% coeficientes_activos_universales$Gen) %>%
  dplyr::select(ModelID, Gen = Gene, metilacion, TPM) %>%
  mutate(
    metilacion = as.numeric(as.character(metilacion)),
    TPM = as.numeric(as.character(TPM))
  ) %>%
  filter(!is.na(metilacion) & !is.na(TPM)) %>%
  filter(metilacion >= 0.70 & TPM < umbral_matematico_exacto) %>% #Hipermetilados y silenciados
  group_by(Gen) %>%
  summarise(Muestras_Con_Silenciamiento_Puro = n(), .groups = "drop")


farma_universales_filtrado <- coeficientes_activos_universales %>%
  inner_join(perfil_muestras_universales_reales, by = "Gen") %>%
  filter(Muestras_Con_Silenciamiento_Puro > 0)


# Exclusión de mecanismos fisiológicos generales o no oncológicos
moa_excluir_univ <- c("ADRENERGIC", "SEROTONIN", "DIURETIC", "NEUROMUSCULAR", "BLOOD PRESSURE", 
                      "PSYCHOTIC", "ANESTHETIC", "OPIOID", "CELL WALL", "BACTERIA", "BACTERIAL", 
                      "ANTIFUNGAL", "ANTIOXIDANT", "ANTIVIRAL", "CALCIUM CHANNEL", "GLUTAMATE", 
                      "HISTAMINE", "DOPAMINE", "ACETYLCHOLINE", "SECRETAC", "GROWTH HORMONE", 
                      "NEPRILYSIN", "CANNABINOID", "SENSING", "UROTENSIN", "PROGESTERONE", 
                      "CCK", "ANDROGEN")

farma_universales_filtrado <- farma_universales_filtrado %>%
  filter(!str_detect(MOA, regex(str_c(moa_excluir_univ, collapse = "|"), ignore_case = TRUE)))


# Se aíslan los 5 MOA más fuertes
mecanismos_ganadores_universales <- farma_universales_filtrado %>% 
  group_by(MOA) %>% 
  summarise(Fuerza_Efecto_Media = mean(abs(Coeficiente), na.rm = TRUE), .groups = "drop") %>% 
  arrange(desc(Fuerza_Efecto_Media), MOA) %>% 
  slice_head(n = 5) %>% 
  pull(MOA)


tabla_frecuencias_universales <- farma_universales_filtrado %>% 
  filter(MOA %in% mecanismos_ganadores_universales) %>%
  mutate(Efecto = if_else(Coeficiente < 0, "Metilación induce sensibilidad", "Metilación induce resistencia")) %>%
  group_by(MOA, Efecto) %>% 
  summarise(Genes_Unicos = n_distinct(Gen), Farmacos_Unicos = n_distinct(Farmaco), 
            Metrica_Normalizada = Genes_Unicos / Farmacos_Unicos, .groups = "drop")


# Traduccion
objeto_grafico_universal_fijo <- tabla_frecuencias_universales %>%
  mutate(MOA_Esp = case_when(
    str_detect(MOA, regex("HISTONE CHAPERONE", ignore_case = TRUE))           ~ "Inhibidor de chaperonas de histonas",
    str_detect(MOA, regex("ESTROGEN RECEPTOR DEGRADER", ignore_case = TRUE))  ~ "Degradador del receptor de estrógenos (SERD)",
    str_detect(MOA, regex("HSP ANTAGONIST|HSP INHIBITOR", ignore_case = TRUE)) ~ "Inhibidor / Antagonista de HSP",
    str_detect(MOA, regex("AURORA KINASE", ignore_case = TRUE))               ~ "Inhibidor de Aurora quinasa",
    str_detect(MOA, regex("UBIQUITIN ACTIVATING ENZYME", ignore_case = TRUE)) ~ "Inhibidor de la enzima activadora de ubiquitina",
    TRUE ~ as.character(MOA) 
  )) %>%
  mutate(Frecuencia_Grafico = if_else(Efecto == "Metilación induce resistencia", Metrica_Normalizada, -Metrica_Normalizada)) %>%
  group_by(MOA_Esp) %>% 
  mutate(Total_Orden = sum(abs(Frecuencia_Grafico))) %>% 
  ungroup()


# Visualización
grafico_moa_universales_espejo_definitivo <- ggplot(objeto_grafico_universal_fijo, 
                                                    aes(x = reorder(MOA_Esp, Total_Orden), 
                                                        y = Frecuencia_Grafico, fill = Efecto)) +
  geom_bar(stat = "identity", alpha = 0.85, width = 0.65) + 
  coord_flip() + 
  scale_fill_manual(values = c("Metilación induce sensibilidad" = "#2B6CB0", 
                               "Metilación induce resistencia"   = "#C53030")) +
  scale_y_continuous(labels = abs, limits = c(-10.0, 10.0)) + 
  theme_minimal(base_size = 14) + 
  labs(
    title = "Mecanismos de acción farmacológica globales",
    x = "Mecanismo de acción clínico (MOA)",
    y = "Densidad de interacciones (Genes predictores por fármaco único)",
    fill = NULL
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 19, color = "#1A202C"),
    axis.title = element_text(face = "bold", size = 19),
    axis.text.y = element_text(size = 17, color = "#2D3748"),
    axis.text.x = element_text(size = 17, color = "#2D3748"), 
    legend.text = element_text(size = 17), 
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

print(grafico_moa_universales_espejo_definitivo)


# 15.1. REPORTE 
# --------------

# Tradución
farma_universales_traducido <- farma_universales_filtrado %>%
  filter(MOA %in% mecanismos_ganadores_universales) %>%
  mutate(Mecanismo_Grafico = case_when(
    str_detect(MOA, regex("HISTONE CHAPERONE", ignore_case = TRUE))           ~ "Inhibidor de chaperonas de histonas",
    str_detect(MOA, regex("ESTROGEN RECEPTOR DEGRADER", ignore_case = TRUE))  ~ "Degradador del receptor de estrógenos (SERD)",
    str_detect(MOA, regex("HSP ANTAGONIST|HSP INHIBITOR", ignore_case = TRUE)) ~ "Inhibidor / Antagonista de HSP",
    str_detect(MOA, regex("AURORA KINASE", ignore_case = TRUE))               ~ "Inhibidor de Aurora quinasa",
    str_detect(MOA, regex("UBIQUITIN ACTIVATING ENZYME", ignore_case = TRUE)) ~ "Inhibidor de la enzima activadora de ubiquitina",
    TRUE ~ as.character(MOA) 
  ))


referencia_grafica_universal <- objeto_grafico_universal_fijo %>%
  dplyr::select(Mecanismo_Grafico = MOA_Esp, Efecto_Grafico = Efecto) %>%
  distinct()


# Extracción de metadatos ómicos para las muestras que pasan el filtro
metadatos_omicos_universales <- tabla_modelos %>%
  filter(Gene %in% farma_universales_traducido$Gen) %>%
  dplyr::select(ModelID, Gen = Gene, metilacion, TPM) %>%
  mutate(
    metilacion = as.numeric(as.character(metilacion)),
    TPM = as.numeric(as.character(TPM))
  ) %>%
  filter(!is.na(metilacion) & !is.na(TPM)) %>%
  filter(metilacion >= 0.70 & TPM < umbral_matematico_exacto) %>% #Hipermetilación y silenciamiento
  group_by(Gen) %>%
  summarise( #Calculo de métricas
    Media_Metilacion = mean(metilacion, na.rm = TRUE),
    Media_Expresion_TPM = mean(TPM, na.rm = TRUE),
    Mediana_Metilacion = median(metilacion, na.rm = TRUE),
    Mediana_Expresion_TPM = median(TPM, na.rm = TRUE),
    Lineas_Celulares_Silenciadas = n(), 
    .groups = "drop"
  )


reporte_tfm_universales_final <- farma_universales_traducido %>%
  dplyr::select(Nombre_Gen = Gen, Farmaco, Coeficiente, MOA) %>%
  mutate(Mecanismo_Grafico = case_when(
    str_detect(MOA, regex("HISTONE CHAPERONE", ignore_case = TRUE))           ~ "Inhibidor de chaperonas de histonas",
    str_detect(MOA, regex("ESTROGEN RECEPTOR DEGRADER", ignore_case = TRUE))  ~ "Degradador del receptor de estrógenos (SERD)",
    str_detect(MOA, regex("HSP ANTAGONIST|HSP INHIBITOR", ignore_case = TRUE)) ~ "Inhibidor / Antagonista de HSP",
    str_detect(MOA, regex("AURORA KINASE", ignore_case = TRUE))               ~ "Inhibidor de Aurora quinasa",
    str_detect(MOA, regex("JAK INHIBITOR", ignore_case = TRUE))                ~ "Inhibidor de JAK",
    str_detect(MOA, regex("DNA BINDING AGENT", ignore_case = TRUE))           ~ "Agonista de unión al ADN (Intercalante)",
    str_detect(MOA, regex("UBIQUITIN ACTIVATING ENZYME", ignore_case = TRUE)) ~ "Inhibidor de la enzima activadora de ubiquitina",
    TRUE ~ as.character(MOA) 
  )) %>%
  # Se asigna el efecto clínico según el signo del modelo
  mutate(Efecto_Real_Signo = if_else(Coeficiente < 0, "Metilación induce sensibilidad", "Metilación induce resistencia")) %>%
  inner_join(
    referencia_grafica_universal, 
    by = c("Mecanismo_Grafico", "Efecto_Real_Signo" = "Efecto_Grafico"),
    relationship = "many-to-many"
  ) %>%
  left_join(metadatos_omicos_universales, by = c("Nombre_Gen" = "Gen")) %>%
  mutate(
    Estadio = "Perfil genómico universal",
    Efecto_Clinico = Efecto_Real_Signo
  ) %>%
  dplyr::select(-Efecto_Real_Signo, -MOA) %>%
  distinct(Mecanismo_Grafico, Nombre_Gen, Farmaco, Efecto_Clinico, .keep_all = TRUE) %>%
  group_by(Mecanismo_Grafico, Efecto_Clinico) %>%
  arrange(Mecanismo_Grafico, Efecto_Clinico, desc(abs(Coeficiente))) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  arrange(desc(abs(Coeficiente))) %>%
  dplyr::select(
    Estadio, Mecanismo_Grafico, Nombre_Gen, Farmaco, Coeficiente, Efecto_Clinico,
    Media_Metilacion, Media_Expresion_TPM, Mediana_Metilacion, Mediana_Expresion_TPM,
    Lineas_Celulares_Silenciadas
  )


write.csv2(reporte_tfm_universales_final, "Reporte_TFM_Universales_farma.csv", row.names = FALSE)

print(as.data.frame(reporte_tfm_universales_final))


# ========================================
# FASE 16: Paquetes y versiones empleadas
# ========================================

sessionInfo()

