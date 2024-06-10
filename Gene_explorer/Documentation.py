import wx
import wx.html2  # Import the wx.html2 module for HTML rendering

def get_cry_gene_information_and_sources():
    """
    Retrieves formatted Cry gene information as HTML, including external links.
    """
    cry_gene_info = """
    <html>
    <body>
    <h2>Understanding Results:</h2>
    <p><strong>E-value:</strong><br>
    Extremely low E-values (approaching 0) and high scores, which are indicative of strong matches.</p>
    
    <p><strong>Score:</strong><br>
    The score in a BLAST alignment is a measure of the similarity between the query sequence and the subject sequence. 
    It is calculated based on the sum of the match and mismatch scores, as well as gap penalties. 
    It is a measure of alignment quality<br>
    <strong>Match score:</strong> Points are added for identical or similar residues nucleotides at each position.<br>
    <strong>Mismatch penalty:</strong> Points are subtracted for non-identical residues.<br>
    <strong>Gap penalties:</strong> Points are subtracted for gaps introduced to optimize the alignment.<br>
    A high score generally indicates a strong alignment between the query and subject sequences, suggesting significant similarity.<br>
    The definition of a "high score" can vary depending on the context (for nucleotide sequences, scores above 350 can be considered high and indicative of significant similarity).</p>
    
    <p><strong>Length of Alignments:</strong><br>
    Length represents the number of nucleotides in the aligned part of the subject sequence.<br>
    The length values (e.g., 1018) will indicate the number of nucleotides in the aligned part of the subject sequence.</p>
    
    <h2>Cry genes Information:</h2>
    <p>Cry genes encode proteins belonging to the Cry toxin family, commonly found in Bacillus thuringiensis.
    These toxins exhibit insecticidal activity against various insect pests, particularly lepidopterans (butterflies and moths).
    The specific insecticidal properties and target insects vary depending on the Cry gene type. In addition to the Cry toxins, 
    some strains of Bt, like Bt israelensis, produce another toxic crystal, named cytolytic protein or Cyt toxin. 
    This Cyt toxin increases the efficiency of Bt in dipteran insects. As of October 2018, approximately 846 cry and cyt genes had been discovered.
    A few of which have used in genetic mofication of crops and these are the ones that have been included in this software. As research goes,
    the database will be expanded to accoomodate them</p>
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
