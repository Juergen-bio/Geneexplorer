import wx
import wx.html2 as webview
from Blastn import run_blast
from Input_output import InputOutputManager
from Documentation import get_cry_gene_information, get_external_sources


class GeneExplorerApp(wx.Frame):
    def __init__(self):
        super().__init__(None, title="Gene Explorer", size=(800, 600))

        panel = wx.Panel(self)

        # E-value Threshold
        evalue_panel = wx.Panel(panel)
        evalue_label = wx.StaticText(evalue_panel, label="E-value:")
        self.evalue_text = wx.TextCtrl(evalue_panel, value="0.001", size=(80, -1))

        evalue_sizer = wx.BoxSizer(wx.HORIZONTAL)
        evalue_sizer.Add(evalue_label, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        evalue_sizer.Add(self.evalue_text, 0, wx.ALL, 5)
        evalue_panel.SetSizer(evalue_sizer)

        # Select Query Sequence
        self.query_seq_label = wx.StaticText(panel, label="Select Query Sequence File:")
        self.query_seq_text = wx.TextCtrl(panel, size=(400, -1), style=wx.TE_READONLY)
        self.select_query_button = wx.Button(panel, label="Select Query Sequence")
        self.select_query_button.Bind(wx.EVT_BUTTON, self.on_select_query)

        # Select Cry Genes Directory
        self.cry_genes_label = wx.StaticText(panel, label="Select Directory containing Cry Gene Files:")
        self.cry_genes_text = wx.TextCtrl(panel, size=(400, -1), style=wx.TE_READONLY)
        self.select_cry_genes_button = wx.Button(panel, label="Select Cry Genes Directory")
        self.select_cry_genes_button.Bind(wx.EVT_BUTTON, self.on_select_cry_genes)

        # Run BLASTN Button
        self.run_blast_button = wx.Button(panel, label="Run BLASTn")
        self.run_blast_button.Bind(wx.EVT_BUTTON, self.on_run_blast)
        self.run_blast_button.Disable()  # Initially disable the button

        # Result Text
        self.result_text = wx.TextCtrl(panel, style=wx.TE_MULTILINE | wx.TE_READONLY, size=(780, 300))
        self.status_label = wx.StaticText(panel, label="Status: Ready")

        # Cry Gene Information Panel
        cry_gene_info = get_cry_gene_information()
        self.cry_gene_info_webview = webview.WebView.New(panel)
        self.cry_gene_info_webview.SetPage(cry_gene_info, "")

        # External Links Panel
        external_links_html = get_external_sources()
        self.external_links_webview = webview.WebView.New(panel)
        self.external_links_webview.SetPage(external_links_html, "")

        # Layout using sizers
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(evalue_panel, 0, wx.ALL | wx.EXPAND, 10)
        sizer.Add(self.query_seq_label, 0, wx.ALL, 10)
        sizer.Add(self.query_seq_text, 0, wx.ALL | wx.EXPAND, 10)
        sizer.Add(self.select_query_button, 0, wx.ALL | wx.CENTER, 10)
        sizer.Add(self.cry_genes_label, 0, wx.ALL, 10)
        sizer.Add(self.cry_genes_text, 0, wx.ALL | wx.EXPAND, 10)
        sizer.Add(self.select_cry_genes_button, 0, wx.ALL | wx.CENTER, 10)
        sizer.Add(self.run_blast_button, 0, wx.ALL | wx.CENTER, 10)
        sizer.Add(self.status_label, 0, wx.ALL, 10)
        sizer.Add(self.result_text, 1, wx.ALL | wx.EXPAND, 10)
        sizer.Add(wx.StaticText(panel, label="Cry Gene Information:"), 0, wx.ALL, 10)
        sizer.Add(self.cry_gene_info_webview, 0, wx.ALL | wx.EXPAND, 10)
        sizer.Add(wx.StaticText(panel, label="External Links:"), 0, wx.ALL, 10)
        sizer.Add(self.external_links_webview, 0, wx.ALL | wx.EXPAND, 10)

        panel.SetSizer(sizer)

        self.selected_query_path = None
        self.selected_cry_genes_path = None

        self.io_manager = InputOutputManager(self)  # Create an instance of InputOutputManager

    def on_select_query(self, event):
        self.io_manager.on_select_query(event)

    def on_select_cry_genes(self, event):
        self.io_manager.on_select_cry_genes(event)

    def on_run_blast(self, event):
        if not self.selected_query_path:
            wx.MessageBox("Please select the query sequence file.", "Error", wx.OK | wx.ICON_ERROR)
            return
        if not self.selected_cry_genes_path:
            wx.MessageBox("Please select the directory containing Cry gene files.", "Error", wx.OK | wx.ICON_ERROR)
            return

        try:
            run_blast(self.selected_query_path, self.selected_cry_genes_path,
                      float(self.evalue_text.GetValue()), self.result_text, self.status_label, self.run_blast_button)

            self.status_label.SetLabel("Status: Completed")

        except Exception as e:
            wx.MessageBox(f"An error occurred: {str(e)}", "Error", wx.OK | wx.ICON_ERROR)
            self.status_label.SetLabel("Status: Error")


if __name__ == "__main__":
    app = wx.App()
    frame = GeneExplorerApp()
    frame.Show()
    app.MainLoop()
