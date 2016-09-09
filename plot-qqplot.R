#################################
#QQplot

#Plots quantile-quantile based on log-10 P-values. 

#Requires:
#Input file with data to read with p-values, at the moment only reads Plink's output from quantitative trait analysis 

#Antonio Berlanga-Taylor
#26/01/2016

#################################

args = commandArgs(trailingOnly = TRUE)

#Set input and output files:
input_file = as.character(args[1]) #'P140343-Results_FinalReport_clean-individuals_and_SNPs'

output_path = file.path(getwd(), paste('qqplot_log10p_values_', input_file, '.png', sep = ''))

# Set column with p-values:
p_value_col = as.numeric(args[2])

# Read data:
data = read.table(input_file, header = TRUE)
#head(data)

# Calculate observed vs expected -log10 p-values:
#obs = -log10(sort(data[data$TEST == 'ADD', ]$P))
data = data[which(!is.na(data[, p_value_col])), ]
data[, p_value_col] = as.numeric(as.character(data[, p_value_col]))

obs = -log10(sort(data[, p_value_col]))

#summary(obs)
exp = -log10(c(1:length(obs)) / (length(obs) + 1 ))
#summary(exp) 

#head(data)
#dim(data)

#Plot:
png(output_path)
plot(exp, obs, ylab = 'Observed (-logP)', xlab = 'Expected(-logP)', ylim = c(0, 20), xlim = c(0, 7))
lines(c(0, 7), c(0, 7), col = 1, lwd = 2)
dev.off()
 

q()
