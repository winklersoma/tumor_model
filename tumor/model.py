import tumor

from mesa.model import Model
from mesa.space import MultiGrid
from mesa.time import RandomActivation
from mesa.datacollection import DataCollector

import random


def metastasis_happening():
    """
    Determines whether metastasis will be happening in the current step.
    :return bool: True if the stem cell will divide symmetrically (into two stem cells)
    """
    value = random.randint(0, 50)
    if value == 0:
        return True


class TumorModel(Model):
    def __init__(self, width: int, height: int):
        super().__init__()
        self.schedule = RandomActivation(model=self)
        self.grid = MultiGrid(width=width, height=height,
                              torus=True)
        a = tumor.TumorCell(unique_id=0,
                            model=self, stem=True, metastatic=False)
        self.schedule.add(agent=a)
        # placing the stem cell in the center of the grid
        x = int(height / 2)
        # x means the vertical axis here (rotated cartesian)
        y = int(width / 2)
        self.grid.place_agent(agent=a,
                              pos=(x, y))
        # unique_id counter for later use
        self.num_id = 1
        # cell counter for visualization
        self.cell_number = 1
        self.datacollector = DataCollector(model_reporters={
            "num_cells": cell_counter
        })
        self.datacollector.collect(model=self)  # saves initial state

    def step(self):
        self.schedule.step()
        for agent in self.schedule.agents:
            agent: tumor.TumorCell  # technical code (unnecessary for running)
            # metastasis (symmetric division of stem cell)
            self.metastasis(agent)
            # asymmetric division
            self.division(agent)
            # death
            self.death(agent)
        self.datacollector.collect(model=self)
        self.cell_number = cell_counter(model=self)

    def death(self, agent):
        """
        Method for removing the old cells from the model.
        :param agent: the cells to whom the death could happen
        :return: None
        """
        if agent.age >= 7 and not (agent.stem or agent.metastatic):
            self.schedule.remove(agent=agent)
            self.grid.remove_agent(agent=agent)

    def division(self, agent):
        """
        Method for the asymmetric division of the cells into transitory cells.
        :param agent: the cells to whom the division could happen
        :return: None
        """
        if agent.stem or agent.metastatic or (agent.age < 3 and agent.divide == 1):
            a = tumor.TumorCell(unique_id=(self.num_id + 1),
                                model=self, stem=False, metastatic=False)
            self.num_id += 1
            self.place_cell(agent=agent, cell=a)

    def metastasis(self, agent):
        """
        Method for the symmetric division of the stem cell creating another moving stem.
        :param agent: the cells to whom the division could happen
        :return: None
        """
        if agent.stem and metastasis_happening():
            a = tumor.TumorCell(unique_id=(self.num_id + 1), model=self, stem=False, metastatic=True)
            self.num_id += 1
            self.place_cell(agent=agent, cell=a)

    def place_cell(self, agent, cell):
        """
        Method for placing the new cell around the existing agent.
        :param agent: The already existing cell around which we place the new
        :param cell: The cell to place
        :return: None
        """
        self.schedule.add(agent=cell)
        cells_to_place = self.grid.get_neighborhood(
            pos=agent.pos,
            moore=False,
            include_center=False,
            radius=1)
        destination_cell = self.random.choice(cells_to_place)
        self.grid.place_agent(agent=cell,
                              pos=destination_cell)


def cell_counter(model: TumorModel):
    """
    Method for counting the number of living tumor cells.
    :param model: The model of which the agents we aim to count.
    :return: the number of agents/cells.
    """
    return model.schedule.get_agent_count()
