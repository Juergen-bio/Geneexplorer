import wx
import wx.html2 as webview
from Blastn import run_blast
from Input_output import InputOutputManager
from Documentation import get_cry_gene_information_and_sources
from result import ResultHandler

class GeneExplorerApp(wx.Frame):
    def __init__(self):
        super().__init__(None, title="Gene Explorer", size=(800, 800))
        panel = wx.Panel(self)
        notebook = wx.Notebook(panel)

        # Query Setup Panel
        query_panel = wx.Panel(notebook)
        self.create_query_setup_ui(query_panel)
        notebook.AddPage(query_panel, "Query Setup")

        # Results Panel
        results_panel = wx.Panel(notebook)
        self.create_results_ui(results_panel)
        notebook.AddPage(results_panel, "Results")

        # Cry Gene Information Panel
        info_panel = wx.Panel(notebook)
        self.create_info_ui(info_panel)
        notebook.AddPage(info_panel, "Documentation")

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(notebook, 1, wx.EXPAND)
        panel.SetSizer(sizer)

        self.selected_query_path = None
        self.selected_cry_genes_path = None
        self.io_manager = InputOutputManager(self)

    def create_query_setup_ui(self, panel):
        # E-value, Word Size, and Gap Penalties Panel
        params_panel = wx.Panel(panel)
        evalue_label = wx.StaticText(params_panel, label="E-value:")
        self.evalue_text = wx.TextCtrl(params_panel, value="0.001", size=(80, -1))
        
        word_size_label = wx.StaticText(params_panel, label="Word Size:")
        self.word_size_text = wx.TextCtrl(params_panel, value="11", size=(80, -1))
        
        gap_penalties_label = wx.StaticText(params_panel, label="Gap Penalties:")
        self.gap_open_text = wx.TextCtrl(params_panel, value="5", size=(40, -1))
        self.gap_extend_text = wx.TextCtrl(params_panel, value="2", size=(40, -1))

        params_sizer = wx.BoxSizer(wx.HORIZONTAL)
        params_sizer.Add(evalue_label, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        params_sizer.Add(self.evalue_text, 0, wx.ALL, 5)
        params_sizer.Add(word_size_label, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        params_sizer.Add(self.word_size_text, 0, wx.ALL, 5)
        params_sizer.Add(gap_penalties_label, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        params_sizer.Add(wx.StaticText(params_panel, label="Open:"), 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        params_sizer.Add(self.gap_open_text, 0, wx.ALL, 5)
        params_sizer.Add(wx.StaticText(params_panel, label="Extend:"), 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        params_sizer.Add(self.gap_extend_text, 0, wx.ALL, 5)
        params_panel.SetSizer(params_sizer)

        self.query_seq_label = wx.StaticText(panel, label="Select Query Sequence File:")
        self.query_seq_text = wx.TextCtrl(panel, size=(400, -1), style=wx.TE_READONLY)
        self.select_query_button = wx.Button(panel, label="Select Query Sequence")
        self.select_query_button.Bind(wx.EVT_BUTTON, self.on_select_query)

        self.cry_genes_label = wx.StaticText(panel, label="Select Directory containing Cry Gene Files:")
        self.cry_genes_text = wx.TextCtrl(panel, size=(400, -1), style=wx.TE_READONLY)
        self.select_cry_genes_button = wx.Button(panel, label="Select Cry Genes Directory")
        self.select_cry_genes_button.Bind(wx.EVT_BUTTON, self.on_select_cry_genes)

        self.run_blast_button = wx.Button(panel, label="Run BLASTn")
        self.run_blast_button.Bind(wx.EVT_BUTTON, self.on_run_blast)
        self.run_blast_button.Disable()

        query_sizer = wx.BoxSizer(wx.VERTICAL)
        query_sizer.Add(params_panel, 0, wx.ALL | wx.EXPAND, 10)
        query_sizer.Add(self.query_seq_label, 0, wx.ALL, 10)
        query_sizer.Add(self.query_seq_text, 0, wx.ALL | wx.EXPAND, 10)
        query_sizer.Add(self.select_query_button, 0, wx.ALL | wx.CENTER, 10)
        query_sizer.Add(self.cry_genes_label, 0, wx.ALL, 10)
        query_sizer.Add(self.cry_genes_text, 0, wx.ALL | wx.EXPAND, 10)
        query_sizer.Add(self.select_cry_genes_button, 0, wx.ALL | wx.CENTER, 10)
        query_sizer.Add(self.run_blast_button, 0, wx.ALL | wx.CENTER, 10)
        panel.SetSizer(query_sizer)

    def create_results_ui(self, panel):
        self.result_text = wx.TextCtrl(panel, style=wx.TE_MULTILINE | wx.TE_READONLY, size=(780, 400))
        self.status_label = wx.StaticText(panel, label="Status: Ready")
        self.progress_bar = wx.Gauge(panel, range=100, size=(780, 25))

        results_sizer = wx.BoxSizer(wx.VERTICAL)
        results_sizer.Add(self.status_label, 0, wx.ALL, 10)
        results_sizer.Add(self.result_text, 1, wx.ALL | wx.EXPAND, 10)
        results_sizer.Add(self.progress_bar, 0, wx.ALL | wx.EXPAND, 10)
        panel.SetSizer(results_sizer)

    def create_info_ui(self, panel):
        cry_gene_info_and_sources = get_cry_gene_information_and_sources()
        self.cry_gene_info_webview = webview.WebView.New(panel, size=(780, 600))
        self.cry_gene_info_webview.SetPage(cry_gene_info_and_sources, "")

        info_sizer = wx.BoxSizer(wx.VERTICAL)
        info_sizer.Add(self.cry_gene_info_webview, 1, wx.ALL | wx.EXPAND, 10)
        panel.SetSizer(info_sizer)

    def on_select_query(self, event):
        self.io_manager.on_select_query(event)

    def on_select_cry_genes(self, event):
        self.io_manager.on_select_cry_genes(event)

    def on_run_blast(self, event):
        print(f"Query Path: {self.selected_query_path}")  # Debug print
        print(f"Cry Genes Path: {self.selected_cry_genes_path}")  # Debug print
        if not self.selected_query_path:
            wx.MessageBox("Please select the query sequence file.", "Error", wx.OK | wx.ICON_ERROR)
            return
        if not self.selected_cry_genes_path:
            wx.MessageBox("Please select the directory containing Cry gene files.", "Error", wx.OK | wx.ICON_ERROR)
            return

        self.progress_bar.Pulse()
        try:
            result_handler = ResultHandler(self.result_text, self.status_label, self.progress_bar, self.run_blast_button)
            run_blast(self.selected_query_path, self.selected_cry_genes_path,
                      float(self.evalue_text.GetValue()), int(self.word_size_text.GetValue()), 
                      int(self.gap_open_text.GetValue()), int(self.gap_extend_text.GetValue()), 
                      result_handler)
        except Exception as e:
            result_handler.handle_error(str(e))

if __name__ == "__main__":
    app = wx.App()
    frame = GeneExplorerApp()
    frame.Show()
    app.MainLoop()
