class Graph:
    def __init__(self, name):
        self.name = name
        self.lines = []
        self.groups = []
        self.traceArray = []

    def getTraceArray(self):
        return f"[{','.join(self.traceArray)}]"

    def addData(self, varName, X, Y, rgbArray = (130,130,130), markerSize = 5):
        name = self.name + '_' + varName
        out = f"""
var {name} = {{
  x: {X},
  y: {Y},
  mode: 'markers',
  type: 'scatter',
  marker: {{ size: {markerSize}, color: 'rgb{rgbArray}' }}

}};
"""
        self.groups.append(out)
        self.traceArray.append(name)
        return out

    def addLine(self, varName, p1, p2, rgbArray = (250, 200, 200), lineWidth = 1):
        name = self.name + '_' + varName
        out = f"""
var {name} = {{
  x: [{p1[0]},{p2[0]}],
  y: [{p1[1]},{p2[1]}],
  mode: 'lines',
  type: 'scatter',
  name: 'nGenes',
  line: {{
    color: 'rgb{rgbArray}',
    width: {lineWidth}
  }}
}};
"""
        self.groups.append(out)
        self.traceArray.append(name)
        return out

    def getGraphData(self):
        return "\n".join(self.groups) + "\n".join(self.lines)

    def addMetadata(self, metadata):
        self.metadata = metadata

    def getMetadata(self):
        return "\t".join([f"{key}={round(value, 2)}" for key,value in self.metadata.items()])
