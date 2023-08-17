from mesa.visualization.modules import CanvasGrid
from mesa.visualization.ModularVisualization import ModularServer

import src.tumor as tumor


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


canvas_element = CanvasGrid(tumor_model_portrayal, 40, 40, 500, 500)

model_params = {
    "height": 40,
    "width": 40,
    "num_agent": 1
}

server = ModularServer(
    tumor.TumorModel, [canvas_element], "Tumor model", model_params
)
