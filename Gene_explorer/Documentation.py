import wx  # Import the wx module for GUI components
import wx.html2  # Import the wx.html2 module for HTML rendering

def get_cry_gene_information():
    """
    Retrieves formatted Cry gene information as HTML.
    """
    cry_gene_info = """
    <html>
    <body>
    <h2>Cry Gene Information:</h2>
    <p>Cry genes encode proteins belonging to the Cry toxin family, commonly found in Bacillus thuringiensis.
    These toxins exhibit insecticidal activity against various insect pests, particularly lepidopterans (butterflies and moths).
    The specific insecticidal properties and target insects vary depending on the Cry gene type.</p>
    </body>
    </html>
    """
    return cry_gene_info


def get_external_sources():
    """
    Retrieves formatted HTML content for external links related to Cry genes.
    """
    external_links = [
        "https://en.wikipedia.org/wiki/Cry_toxin",
        "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6148975/",
        "https://www.sciencedirect.com/topics/agricultural-and-biological-sciences/cry-proteins"
    ]
    
    links_html = "<html><body><h2>External Links:</h2><ul>"
    for link in external_links:
        links_html += f'<li><a href="{link}">{link}</a></li>'
    links_html += "</ul></body></html>"
    
    return links_html

if __name__ == "__main__":
    # Example usage if executed as a standalone script
    app = wx.App()
    frame = wx.Frame(None, title="Documentation Example")

    # Retrieve Cry gene information HTML
    cry_gene_info_html = get_cry_gene_information()

    # Retrieve external links HTML
    external_links_html = get_external_sources()

    # Create wx.html2.WebView for displaying Cry gene information
    cry_gene_webview = wx.html2.WebView.New(frame)
    cry_gene_webview.SetPage(cry_gene_info_html, "")

    # Create wx.html2.WebView for displaying external links
    external_links_webview = wx.html2.WebView.New(frame)
    external_links_webview.SetPage(external_links_html, "")

    # Layout the frame
    sizer = wx.BoxSizer(wx.VERTICAL)
    sizer.Add(cry_gene_webview, 1, wx.EXPAND | wx.ALL, 10)
    sizer.Add(external_links_webview, 1, wx.EXPAND | wx.ALL, 10)
    frame.SetSizer(sizer)

    frame.SetSize((800, 600))
    frame.Show()
    app.MainLoop()
