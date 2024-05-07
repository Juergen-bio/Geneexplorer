import os
import wx
import threading
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO

def run_blast(selected_query_path, selected_cry_genes_path, evalue_threshold, result_text, status_label, run_blast_button):
    try:
        query_sequences = list(SeqIO.parse(selected_query_path, "fasta"))

        if not query_sequences:
            wx.MessageBox("No sequences found in the selected query file.", "Error", wx.OK | wx.ICON_ERROR)
            return

        result_text.SetValue("")  # Clear previous results
        status_label.SetLabel("Status: Running BLASTn...")
        run_blast_button.Disable()  # Disable button during BLAST

        def process_blast_results(query_sequences, cry_genes_path, evalue_threshold, result_text):
            for file_name in os.listdir(cry_genes_path):
                if file_name.endswith(".fasta"):
                    cry_gene_path = os.path.join(cry_genes_path, file_name)
                    cry_gene_name = os.path.basename(cry_gene_path).replace(".fasta", "")
                    cry_gene_content = str(SeqIO.read(cry_gene_path, "fasta").seq)

                    for seq_record in query_sequences:
                        query_gene = str(seq_record.seq)

                        result_handle = NCBIWWW.qblast("blastn", "nr", query_gene, expect=evalue_threshold)
                        blast_records = NCBIXML.parse(result_handle)

                        cry_gene_found = False

                        for record in blast_records:
                            for alignment in record.alignments:
                                result_text.AppendText(f"Query Sequence: {seq_record.id}\n")
                                result_text.AppendText(f"Alignment Title: {alignment.title}\n")
                                result_text.AppendText(f"Length: {alignment.length}\n")
                                result_text.AppendText(f"E-value: {alignment.hsps[0].expect}\n")
                                result_text.AppendText(f"Score: {alignment.hsps[0].score}\n\n")

                                # Check if Cry gene is detected
                                if cry_gene_content in alignment.title:
                                    cry_gene_found = True

                        # Provide specific message based on Cry gene detection
                        if cry_gene_found:
                            result_text.AppendText(f"Plant is likely a GMO (contains {cry_gene_name})\n")
                        else:
                            result_text.AppendText(f"Plant doesn't contain {cry_gene_name}\n")

            status_label.SetLabel("Status: Completed")
            run_blast_button.Enable()  # Re-enable button after BLAST

        blast_thread = threading.Thread(target=process_blast_results,
                                         args=(query_sequences, selected_cry_genes_path, evalue_threshold, result_text))
        blast_thread.start()

    except Exception as e:
        wx.MessageBox(f"An error occurred: {str(e)}", "Error", wx.OK | wx.ICON_ERROR)
        status_label.SetLabel("Status: Error")
