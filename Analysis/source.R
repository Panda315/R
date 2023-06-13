totalDark24_data <- read.csv("/Users/panda/Downloads/dark24real.csv", header = TRUE, sep = ",", row.names = 1)
totalDark24_data <- data.frame(totalDark24_data)
totalLight24_data <- read.csv("/Users/panda/Downloads/light24.csv", header = TRUE, sep = ",", row.names = 1)
totalLight_data <- data.frame(totalLight24_data)
totalLightDark12_data <- read.csv("/Users/panda/Downloads/lightdark1212.csv", header = TRUE, sep = ",", row.names = 1)
totalLightDark12_data <- data.frame(totalLightDark12_data)
meta_data <- read.csv("/Users/panda/Downloads/24-12-0-CD-ABX.csv", header = TRUE, sep = ",", row.names = 1)
meta_data <- data.frame(meta_data)

totalDark24_data <- t(totalDark24_data)
dds <- DESeqDataSetFromMatrix(countData = totalDark24_data, colData = totalLightDark12_data)
dds
