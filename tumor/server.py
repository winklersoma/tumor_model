import tumor as tumor

from mesa.visualization.ModularVisualization import ModularServer
from mesa.visualization.modules import CanvasGrid, TextElement, ChartModule


class CounterElement(TextElement):
    """
    Display a text count of how many cells there are.
    """
    def render(self, model):
        return "Cancer cells: " + str(model.cell_number)


def tumor_model_portrayal(person):
    if person is None:
        return

    portrayal = {}
    if person.stem or person.metastatic:
        portrayal["Shape"] = "pics/stem_blue.png"
        portrayal["scale"] = 0.9
        portrayal["Layer"] = 4

    elif (person.stem and person.metastatic) is False and person.age < 3:
        portrayal["Shape"] = "pics/transitory_red.png"
        portrayal["scale"] = 0.9
        portrayal["Layer"] = 3

    elif (person.stem and person.metastatic) is False and 4 <= person.age < 6:
        portrayal["Shape"] = "pics/transitory_white.png"
        portrayal["scale"] = 0.9
        portrayal["Layer"] = 1

    elif (person.stem and person.metastatic) is False and 6 <= person.age:
        portrayal["Shape"] = "pics/transitory_black.png"
        portrayal["scale"] = 0.9
        portrayal["Layer"] = 2

    else:
        pass

    return portrayal


cell_chart = ChartModule([{"Label": "num_cells", "Color": "Black"}])
counter_element = CounterElement()
canvas_element = CanvasGrid(tumor_model_portrayal, 40, 40, 500, 500)

model_params = {"height": 40, "width": 40}

server = ModularServer(tumor.TumorModel, [canvas_element, counter_element, cell_chart], "Tumor model", model_params)
