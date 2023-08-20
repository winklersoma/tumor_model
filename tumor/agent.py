import tumor as tumor

from mesa.agent import Agent

import random


def division() -> bool:
    """
    Decides if the cell is to divide in the current step or not.
    :return bool: True if the cell is to divide.
    """
    if random.randint(0, 2) == 0:
        return True


class TumorCell(Agent):
    """
    Constructor for the TumorCell Agent.
    """
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
        self.movement()

    def movement(self):
        """
        Movement method for the tumor cells.
        :return: None
        """
        self.model: tumor.TumorModel
        # using the moore neighbourhood to get destination options
        cells_to_move = self.model.grid.get_neighborhood(
            pos=self.pos,
            moore=True,
            include_center=False,
            radius=1)
        # movement for transitory cells
        if (not self.stem) and (not self.metastatic):
            destination_cell = self.model.random.choice(cells_to_move)
        # movement for the young metastatic stem cell
        if self.metastatic and self.age < (min(self.model.grid.width, self.model.grid.height) / 3):
            destination_cell = cells_to_move[self.unique_id % 8]  # The direction of the movement is constant!
        # original stem and old metastatic stems stay put
        else:
            destination_cell = self.pos
        self.model.grid.move_agent(agent=self, pos=destination_cell)
