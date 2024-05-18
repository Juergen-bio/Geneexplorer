import wx
import wx.html2  # Import the wx.html2 module for HTML rendering

def get_cry_gene_information_and_sources():
    """
    Retrieves formatted Cry gene information as HTML, including external links.
    """
    cry_gene_info = """
    <html>
    <body>
    <h2>Cry Gene Information:</h2>
    <p>Cry genes encode proteins belonging to the Cry toxin family, commonly found in Bacillus thuringiensis.
    These toxins exhibit insecticidal activity against various insect pests, particularly lepidopterans (butterflies and moths).
    The specific insecticidal properties and target insects vary depending on the Cry gene type.</p>
    <h2>External Links:</h2>
    <ul>
        <li><a href="https://en.wikipedia.org/wiki/Cry_toxin">https://en.wikipedia.org/wiki/Cry_toxin</a></li>
        <li><a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6148975/">https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6148975/</a></li>
        <li><a href="https://www.sciencedirect.com/topics/agricultural-and-biological-sciences/cry-proteins">https://www.sciencedirect.com/topics/agricultural-and-biological-sciences/cry-proteins</a></li>
    </ul>
    </body>
    </html>
    """
    return cry_gene_info

if __name__ == "__main__":
    app = wx.App()
    frame = wx.Frame(None, title="Documentation Example")
    cry_gene_info_and_sources_html = get_cry_gene_information_and_sources()
    cry_gene_info_webview = wx.html2.WebView.New(frame)
    cry_gene_info_webview.SetPage(cry_gene_info_and_sources_html, "")

    sizer = wx.BoxSizer(wx.VERTICAL)
    sizer.Add(cry_gene_info_webview, 1, wx.EXPAND | wx.ALL, 10)
    frame.SetSizer(sizer)

    frame.SetSize((800, 800))
    frame.Show()
    app.MainLoop()
