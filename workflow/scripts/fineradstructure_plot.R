#!/usr/bin/env Rscript
# Plot fineRADstructure results
# Visualizes the co-ancestry tree from fineRADstructure MCMC tree output

library(ape)
library(ggtree)
library(ggplot2)
library(treeio)

# Try to load XML package, but handle gracefully if not available
if (!require("XML", quietly = TRUE)) {
  cat("WARNING: XML package not found. Attempting alternative methods...\n")
}

# Prevent creation of Rplots.pdf
pdf(NULL)

# Redirect all output to log file
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

# Get Snakemake inputs
mcmcTree_xml <- snakemake@input[["mcmcTree"]]
chunks_file <- snakemake@input[["chunks"]]

# Get Snakemake outputs
output_pdf <- snakemake@output[["pdf"]]
output_rds <- snakemake@output[["rds"]]

cat("=== fineRADstructure Plotting ===\n")
cat("MCMC Tree XML:", mcmcTree_xml, "\n")
cat("Chunks file:", chunks_file, "\n")
cat("Output PDF:", output_pdf, "\n")
cat("Output RDS:", output_rds, "\n")
cat("===============================\n\n")

# Function to extract tree from fineRADstructure XML
# fineRADstructure uses fineSTRUCTURE's own XML format: the tree is stored as
# Newick text inside the last <outputFile> node (not phyloXML or BEAST).
extract_tree_from_xml <- function(xml_file) {
  cat("Parsing XML file:", xml_file, "\n")

  # Method 1: fineSTRUCTURE native format (official fineRADstructure layout)
  # Tree is Newick in the last outputFile node; see FinestructureLibrary.R extractTree()
  if (require("XML", quietly = TRUE)) {
    tryCatch({
      cat("Attempting to read as fineSTRUCTURE XML (last outputFile = Newick)...\n")
      xml_doc <- xmlParse(xml_file)
      xml_root <- xmlRoot(xml_doc)
      output_file_nodes <- getNodeSet(xml_root, "//outputFile")
      if (length(output_file_nodes) > 0) {
        last_node <- output_file_nodes[[length(output_file_nodes)]]
        newick_string <- xmlValue(last_node)
        newick_string <- trimws(newick_string)
        if (nchar(newick_string) > 0 && grepl("\\(", newick_string)) {
          cat("Successfully read tree from fineSTRUCTURE outputFile (Newick)\n")
          tree <- read.tree(text = newick_string)
          return(tree)
        }
      }
    }, error = function(e) {
      cat("fineSTRUCTURE outputFile extraction failed:", conditionMessage(e), "\n")
    })
  }

  # Method 2: Generic XML fallback – look for //tree or other Newick-containing nodes
  if (require("XML", quietly = TRUE)) {
    tryCatch({
      cat("Attempting to parse XML and extract Newick from //tree or similar...\n")
      xml_doc <- xmlParse(xml_file)
      xml_root <- xmlRoot(xml_doc)
      tree_nodes <- getNodeSet(xml_root, "//tree")
      if (length(tree_nodes) > 0) {
        for (node in tree_nodes) {
          tree_string <- trimws(xmlValue(node))
          if (nchar(tree_string) > 0 && grepl("\\(", tree_string)) {
            cat("Found Newick tree string in XML\n")
            tree <- read.tree(text = tree_string)
            return(tree)
          }
        }
      }
      tree_attrs <- xpathApply(xml_root, "//tree", xmlAttrs)
      if (length(tree_attrs) > 0) {
        for (attr in tree_attrs) {
          if ("newick" %in% names(attr)) {
            cat("Found Newick tree in XML attributes\n")
            tree <- read.tree(text = attr[["newick"]])
            return(tree)
          }
        }
      }
    }, error = function(e) {
      cat("XML fallback parsing failed:", conditionMessage(e), "\n")
    })
  }
  
  # Method 3: Try to extract Newick using shell command (if finestructure tools available)
  cat("Attempting to extract tree using shell command...\n")
  tryCatch({
    # fineRADstructure may have command-line tools to extract tree
    # Try using grep/sed to extract Newick format from XML
    temp_newick <- tempfile(fileext = ".nwk")
    cmd <- paste("grep -o '<tree[^>]*>.*</tree>'", xml_file, "| sed 's/<[^>]*>//g' | head -1 >", temp_newick)
    system(cmd, ignore.stderr = TRUE)
    
    if (file.exists(temp_newick) && file.info(temp_newick)$size > 0) {
      tree_string <- readLines(temp_newick, n = 1, warn = FALSE)
      if (nchar(tree_string) > 0 && grepl("\\(", tree_string)) {
        cat("Extracted Newick string using shell command\n")
        tree <- read.tree(text = tree_string)
        unlink(temp_newick)
        return(tree)
      }
    }
    unlink(temp_newick)
  }, error = function(e) {
    cat("Shell extraction failed:", conditionMessage(e), "\n")
  })
  
  # Method 4: Try reading file as text and searching for Newick pattern
  cat("Attempting to extract Newick from file text...\n")
  tryCatch({
    file_lines <- readLines(xml_file, warn = FALSE)
    # Look for lines that look like Newick format
    for (line in file_lines) {
      if (grepl("^[^<]*\\(.*\\)[^<]*$", line) && grepl(":", line)) {
        # Extract potential Newick string (remove XML tags)
        newick_candidate <- gsub("<[^>]+>", "", line)
        newick_candidate <- trimws(newick_candidate)
        if (nchar(newick_candidate) > 10 && grepl("\\(", newick_candidate)) {
          cat("Found potential Newick string in file text\n")
          tree <- read.tree(text = newick_candidate)
          return(tree)
        }
      }
    }
  }, error = function(e) {
    cat("Text extraction failed:", conditionMessage(e), "\n")
  })
  
  stop("Could not extract tree from XML file. Tried multiple methods. Please check the XML structure or ensure fineRADstructure tools are available.")
}

# Function to create tree plot
create_tree_plot <- function(tree) {
  cat("Creating tree plot...\n")
  
  # Count number of tips to determine plot size
  n_tips <- length(tree$tip.label)
  cat("Number of tips:", n_tips, "\n")
  
  # Calculate appropriate plot dimensions
  plot_height <- max(6, min(50, n_tips * 0.3))
  plot_width <- max(8, min(20, plot_height * 0.8))
  
  # Create the base plot
  p <- ggtree(tree, layout = "rectangular") +
    geom_tiplab(size = 2.5, hjust = -0.05) +
    theme_tree2()
  
  # Add scale bar
  p <- p + geom_treescale(x = 0, y = 0, width = NULL, offset = 1)
  
  # Adjust x-axis limits to accommodate tip labels
  max_x <- max(p$data$x, na.rm = TRUE)
  p <- p + xlim(NA, max_x * 1.3)
  
  return(list(plot = p, width = plot_width, height = plot_height))
}

# Main plotting workflow
tryCatch({
  # Extract tree from XML
  tree <- extract_tree_from_xml(mcmcTree_xml)
  
  # Create plot
  plot_result <- create_tree_plot(tree)
  p <- plot_result$plot
  plot_width <- plot_result$width
  plot_height <- plot_result$height
  
  # Save as PDF
  cat("Saving PDF plot...\n")
  ggsave(
    output_pdf,
    plot = p,
    width = plot_width,
    height = plot_height,
    units = "in",
    limitsize = FALSE
  )
  
  # Save as RDS
  cat("Saving RDS object...\n")
  saveRDS(p, output_rds)
  
  cat("\n=== fineRADstructure Plotting Complete ===\n")
  cat("Plot dimensions:", plot_width, "x", plot_height, "inches\n")
  cat("Number of tips:", length(tree$tip.label), "\n")
  cat("Output PDF:", output_pdf, "\n")
  cat("Output RDS:", output_rds, "\n")
  cat("==========================================\n")
  
}, error = function(e) {
  cat("ERROR in fineRADstructure plotting:", conditionMessage(e), "\n")
  cat("Creating placeholder plot...\n")
  
  # Create placeholder PDF
  pdf(output_pdf, width = 10, height = 8)
  plot.new()
  text(0.5, 0.7, "fineRADstructure Tree Plot", cex = 1.8, font = 2)
  text(0.5, 0.5, paste("MCMC Tree XML:", mcmcTree_xml), cex = 1.2)
  text(0.5, 0.4, paste("Error:", conditionMessage(e)), cex = 1, col = "red")
  text(0.5, 0.2, paste("Generated:", Sys.time()), cex = 0.9, col = "gray")
  dev.off()
  
  # Create placeholder RDS
  p_placeholder <- ggplot() + 
    annotate("text", x = 0.5, y = 0.5, 
             label = paste("Plot generation failed:", conditionMessage(e)),
             size = 5, color = "red")
  saveRDS(p_placeholder, output_rds)
  
  stop("fineRADstructure plotting failed: ", conditionMessage(e))
})

# Close log file
sink()
sink(type = "message")
close(log_file)

