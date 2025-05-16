# Load necessary libraries
library(dplyr)
library(stringr)
library(paletteer)
library(ggplot2)

# Define the directory containing Bismark reports
report_dir <- "~/methylation/sea_lion/bismark_output/alignment"
report_files <- list.files(report_dir, pattern = "_bismark_bt2_PE_report.txt$", full.names = TRUE)

# Initialize an empty list to store extracted data
report_list <- list()

# Loop through all report files and extract relevant values
for (file in report_files) {
  # Read the file as plain text
  file_content <- readLines(file)
  
  # Extract relevant alignment statistics
  TotalPairs <- as.numeric(str_extract(file_content[grep("Sequence pairs analysed in total:", file_content)], "\\d+"))
  UniqueAlignments <- as.numeric(str_extract(file_content[grep("Number of paired-end alignments with a unique best hit:", file_content)], "\\d+"))
  MappingEfficiency <- as.numeric(str_extract(file_content[grep("Mapping efficiency:", file_content)], "\\d+\\.\\d+"))
  NoAlignment <- as.numeric(str_extract(file_content[grep("Sequence pairs with no alignments under any condition:", file_content)], "\\d+"))
  NotUnique <- as.numeric(str_extract(file_content[grep("Sequence pairs did not map uniquely:", file_content)], "\\d+"))
  Discarded <- as.numeric(str_extract(file_content[grep("Sequence pairs which were discarded because genomic sequence could not be extracted:", file_content)], "\\d+"))
  TotalNumberOfCsAnalysed <- as.numeric(str_extract(file_content[grep("Total number of C's analysed", file_content)], "\\d+"))
  TotalMethylatedCsInCpGCotext <- as.numeric(str_extract(file_content[grep("Total methylated C's in CpG context:", file_content)], "\\d+"))
  TotalMethylatedCsInCHGCotext <- as.numeric(str_extract(file_content[grep("Total methylated C's in CHG context:", file_content)], "\\d+"))
  TotalMethylatedCsInCHHCotext <- as.numeric(str_extract(file_content[grep("Total methylated C's in CHH context:", file_content)], "\\d+"))
  TotalMethylatedCsInUnknownCotext <- as.numeric(str_extract(file_content[grep("Total methylated C's in Unknown context", file_content)], "\\d+"))
  TotalUnmethylatedCsInCpGCotext <- as.numeric(str_extract(file_content[grep("Total unmethylated C's in CpG context", file_content)], "\\d+"))
  TotalUnmethylatedCsInCHGCotext <- as.numeric(str_extract(file_content[grep("Total unmethylated C's in CHG context", file_content)], "\\d+"))
  TotalUnmethylatedCsInCHHCotext <- as.numeric(str_extract(file_content[grep("Total unmethylated C's in CHH context", file_content)], "\\d+"))
  TotalUnmethylatedCsInUnknownCotext <- as.numeric(str_extract(file_content[grep("Total unmethylated C's in Unknown context", file_content)], "\\d+"))
  CmethylatedInCpGcontext <- as.numeric(str_extract(file_content[grep("C methylated in CpG context", file_content)], "\\d+"))
  CmethylatedInCHGcontext <- as.numeric(str_extract(file_content[grep("C methylated in CHG context", file_content)], "\\d+"))
  CmethylatedInCHHcontext <- as.numeric(str_extract(file_content[grep("C methylated in CHH context", file_content)], "\\d+"))
  
  
  # Get sample name
  sample_name <- gsub("_bismark_bt2_PE_report.txt", "", basename(file))
  
  # Store extracted data in a list
  report_list[[sample_name]] <- data.frame(
    Sample = sample_name, 
    TotalPairs = TotalPairs, 
    UniqueAlignments = UniqueAlignments, 
    MappingEfficiency = MappingEfficiency, 
    NoAlignment = NoAlignment, 
    NotUnique = Discarded, 
    Discarded = Discarded,
    TotalNumberOfCsAnalysed=TotalNumberOfCsAnalysed,
    TotalMethylatedCsInCpGCotext=TotalMethylatedCsInCpGCotext,
    TotalMethylatedCsInCHGCotext=TotalMethylatedCsInCHGCotext,
    TotalMethylatedCsInCHHCotext=TotalMethylatedCsInCHHCotext,
    TotalMethylatedCsInUnknownCotext=TotalMethylatedCsInUnknownCotext,
    TotalUnmethylatedCsInCpGCotext=TotalUnmethylatedCsInCpGCotext,
    TotalUnmethylatedCsInCHGCotext=TotalUnmethylatedCsInCHGCotext,
    TotalUnmethylatedCsInCHHCotext=TotalUnmethylatedCsInCHHCotext,
    TotalUnmethylatedCsInUnknownCotext=TotalUnmethylatedCsInUnknownCotext,
    CmethylatedInCpGcontext=CmethylatedInCpGcontext,
    CmethylatedInCHGcontext=CmethylatedInCHGcontext,
    CmethylatedInCHHcontext=CmethylatedInCHHcontext)
}

# Combine all extracted data into a single dataframe
final_dataframe <- bind_rows(report_list)

# Converting methylated C into percentages 
final_dataframe$CmethylatedInCpGcontext <- final_dataframe$CmethylatedInCpGcontext/100
final_dataframe$CmethylatedInCHGcontext <- final_dataframe$CmethylatedInCHGcontext/100
final_dataframe$CmethylatedInCHHcontext <- final_dataframe$CmethylatedInCHHcontext/100

# Save the combined dataframe as a CSV file
write.csv(final_dataframe, "~/methylation/sealion_clock/output_csv/combined_bismark_alignment_summary.csv", row.names = FALSE)

# Create summary statistics
summary_stats <- final_dataframe %>%
  summarise(
    Sample  = n(),
    Mean_Mapping_Efficiency = mean(MappingEfficiency , na.rm = TRUE),
    SD_Mapping_Efficiency = sd(MappingEfficiency , na.rm = TRUE),
    Mean_CpG_Methylation = mean(CmethylatedInCpGcontext, na.rm = TRUE),
    SD_CpG_Methylation = sd(CmethylatedInCpGcontext, na.rm = TRUE),
    Mean_Total_Reads = mean(TotalPairs , na.rm = TRUE),
    SD_Total_Reads = sd(TotalPairs , na.rm = TRUE))

# Print summary statistics
print(summary_stats)

### Create visualizations if there's data

# Mapping efficiency plot
mapping_plot <- ggplot(final_dataframe, aes(x = Sample, y = MappingEfficiency)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Mapping Efficiency Across Samples",
       y = "Mapping Efficiency (%)")

ggsave(plot = mapping_plot, filename = "~/methylation/sealion_clock/output_plot/mapping_plot.png" , width = 6, height = 4, dpi = 300)

name(final_dataframe)

# Methylation context plot
methylation_long <- final_dataframe %>%
    pivot_longer(cols = c(CmethylatedInCHHcontext, CmethylatedInCHGcontext, CmethylatedInCpGcontext),
                 names_to = "Context",
                 values_to = "Methylation")

methylation_plot <- ggplot(methylation_long, aes(x = Sample, y = Methylation, fill = Context)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_paletteer_d("nationalparkcolors::Acadia") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
    labs(title = "Methylation Levels by Context", y = "Methylation (%)")

ggsave(plot = methylation_plot, filename = "~/methylation/sealion_clock/output_plot/methylation_plot.png", dpi = 300)











