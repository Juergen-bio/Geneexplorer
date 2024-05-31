import wx
import matplotlib.pyplot as plt
from io import BytesIO

class ResultHandler:
    def __init__(self, result_text, status_label, progress_bar, run_blast_button):
        self.result_text = result_text
        self.status_label = status_label
        self.progress_bar = progress_bar
        self.run_blast_button = run_blast_button

    def display_results(self, results, evalues):
        self.result_text.SetValue(results)
        self.status_label.SetLabel("Status: Completed")
        self.progress_bar.SetValue(100)
        self.run_blast_button.Enable()

        if evalues:
            hist_bitmap = self.show_histogram(evalues)
            self.result_text.AppendText("\nHistogram of E-values:\n")
            self.result_text.WriteImage(hist_bitmap)

    def handle_error(self, error_message):
        wx.MessageBox(f"An error occurred: {error_message}", "Error", wx.OK | wx.ICON_ERROR)
        self.status_label.SetLabel("Status: Error")
        self.progress_bar.SetValue(0)
        self.run_blast_button.Enable()

    def show_histogram(self, evalues):
        fig, ax = plt.subplots()
        ax.hist(evalues, bins=30, edgecolor='black')
        ax.set_title('E-value Distribution')
        ax.set_xlabel('E-value')
        ax.set_ylabel('Frequency')

        buf = BytesIO()
        fig.savefig(buf, format='png')
        buf.seek(0)
        img = wx.Image(buf, wx.BITMAP_TYPE_PNG)
        return img.ConvertToBitmap()