GMO Status Assessment of Maize and Soybean Samples in Kenya


Project Overview
This project focuses on the bioinformatics-driven assessment of Genetically Modified Organism (GMO) status in maize and soybean samples collected from Kenya.
The primary objective was to accurately classify samples as GMO or non-GMO by detecting the presence of specific Cry genes, which are commonly engineered into these crops for pest resistance. 
This work involved comprehensive genomic data analysis and the development of a custom bioinformatics pipeline.

Problem Statement
The accurate and efficient detection of GMOs is critical for food safety, regulatory compliance, and consumer confidence. Traditional methods can be time-consuming and labor-intensive. 
This project addresses the need for a robust, bioinformatics-based method for biomarker detection and classification of GMO status in agricultural samples.

Methodology
Tools & Technologies Used
•	Programming Language: Python 
o	Key Libraries: [Biopython, Pandas, NumPy]
•	Bioinformatics Software: BLAST+ suite
•	GUI Framework: wxPython
•	Data Formats: FASTA
Project Structure
.
├── code/
│   ├── main.py                # Main script for GMO analysis
│   ├── gui.py                 # wxPython GUI implementation
│   ├── analysis_functions.py  # Core bioinformatics analysis functions
│   └── utils.py               # Utility functions (e.g., data parsing)
├── data/
│   ├── reference_cry_genes.fasta # Curated Cry gene sequences
│   ├── sample_genomes/        # Directory for sample genome files
│   └── results/               # Directory for analysis output
├── documentation.docx         # Original project documentation
└── README.md                  # This file

Results & Impact
•	The developed pipeline provides a robust and efficient method for GMO classification based on specific Cry gene presence.
•	By filtering for approved Cry genes, the approach optimizes processing time and focuses on relevant genetic markers.
•	The project demonstrates practical application of bioinformatics principles, genomic data analysis, and Python programming for biomarker detection.

Future Enhancements 
•	Integrate additional gene targets beyond Cry genes.
•	Implement containerization (e.g., Docker) for easier deployment.
•	Develop a web-based interface for broader accessibility.
•	Incorporate more advanced visualization techniques for results.

Author
Juergen George
LinkedIn Profile: www.linkedin.com/in/juergen-george 


