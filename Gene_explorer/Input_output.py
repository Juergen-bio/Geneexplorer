import wx

class InputOutputManager:
    def __init__(self, parent_frame):
        self.parent_frame = parent_frame
        self.selected_query_path = None
        self.selected_cry_genes_path = None

    def on_select_query(self, event):
        wildcard = "Sequence files (*.fasta)|*.fasta|All files (*.*)|*.*"
        dialog = wx.FileDialog(self.parent_frame, "Select Query Sequence File", wildcard=wildcard, style=wx.FD_OPEN)

        if dialog.ShowModal() == wx.ID_OK:
            self.selected_query_path = dialog.GetPath()
            self.parent_frame.query_seq_text.SetValue(self.selected_query_path)
            self.check_enable_blast_button()  # Check and enable BLASTn button

        dialog.Destroy()

    def on_select_cry_genes(self, event):
        dialog = wx.DirDialog(self.parent_frame, "Select Directory containing Cry Gene Files")

        if dialog.ShowModal() == wx.ID_OK:
            self.selected_cry_genes_path = dialog.GetPath()
            self.parent_frame.cry_genes_text.SetValue(self.selected_cry_genes_path)
            self.check_enable_blast_button()  # Check and enable BLASTn button

        dialog.Destroy()

    def check_enable_blast_button(self):
        #Enable the BLASTn button if both query and Cry genes directory are selected.
        if self.selected_query_path and self.selected_cry_genes_path:
            self.parent_frame.run_blast_button.Enable()
