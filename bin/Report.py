class Report:
    def __init__(self,name = "main_report"):
        self.graphs = []
        self.name = name
    def getHead(self):

        out = f"""<!DOCTYPE html>
<html>
    <head>
        <meta charset="utf-8">
        <title>{self.name}</title>
        <script src="https://cdn.plot.ly/plotly-2.16.1.min.js"></script>
    </head>
    <body>
	    <h1 style="font-family: Sans-serif; color: orange; border: grey;">test</h1>

        """
        for g in self.graphs:
            out += f"\n{g.getMetadata()}\n"
            out += f"""\n<div id="{g.name}" style="width:1000px;height:700px;"></div>\n"""

        out += "\n<script>\n"
        return out
    def getTail(self):
        out  = ""
        for g in self.graphs:
            out += f"\nvar data = {g.getTraceArray()};\nPlotly.newPlot('{g.name}', data);\n"
        out+="\n</script>\n</body>\n</html>\n"
        return out

    def getBody(self):
        return "\n".join(list([g.getGraphData() for g in self.graphs]))

    def write(self, filePath):
        handler = open(filePath, 'w')
        handler.write(f"{self.getHead()}\n{self.getBody()}\n{self.getTail()}")
        handler.close()

    def addGraph(self, graph):
        self.graphs.append(graph)
