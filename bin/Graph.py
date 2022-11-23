class Graph:
    def __init__(self, name):
        self.name = name
        self.lines = []
        self.groups = []
        self.traceArray = []

    def addData(self, varName, X, Y, rgbArray, markerSize):
        out = f"""
var {varName} = {{
  x: {X},
  y: {Y},
  mode: 'markers',
  type: 'scatter'
}};
"""
        self.groups.append(out)
        self.traceArray.append(varName)
        return out

    def addLine(self, varName, p1, p2, rgbArray = (150, 150, 150), lineWidth = 2):

        out = f"""
var {varName} = {{
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
        self.traceArray.append(varName)
        return out
