from mesa.agent import Agent
import tumor as tumor

import random


def division() -> bool:
    """
    Decides if the cell is to divide in the current step or not.
    :return bool: True if the cell is to divide.
    """
    if random.randint(0, 2) == 0:
        return True


class TumorCell(Agent):
    def __init__(self, unique_id: int, model: tumor.TumorModel, stem: bool, metastatic: bool):
        super().__init__(unique_id=unique_id,
                         model=model)
        self.age = 0
        # True if the cell is a stem cell
        self.stem = stem
        # If true, then the cell divides in the current round
        self.divide = False
        # True if the cell is a metastatic (stem) cell
        self.metastatic = metastatic

    def step(self):
        self.model: tumor.TumorModel
        # ageing
        self.age += 1
        # mitosis
        self.divide = division()
        # movement
        if not self.stem and not self.metastatic:
            self.move()
        if self.metastatic and self.age < (min(self.model.grid.width, self.model.grid.height) / 3):
            self.metastatic_movement()

    def move(self):
        self.model: tumor.TumorModel
        cells_to_move = self.model.grid.get_neighborhood(
            pos=self.pos,
            moore=True,
            include_center=False,
            radius=1)
        destination_cell = self.model.random.choice(cells_to_move)
        self.model.grid.move_agent(agent=self, pos=destination_cell)

    def metastatic_movement(self):
        self.model: tumor.TumorModel
        cells_to_move = self.model.grid.get_neighborhood(
            pos=self.pos,
            moore=True,
            include_center=False,
            radius=1)
        destination_cell = cells_to_move[self.unique_id % 8]  # Meaning that the direction of the movement is constant!
        self.model.grid.move_agent(agent=self, pos=destination_cell)
