import src.tumor as tumor

from mesa.model import Model
from mesa.space import MultiGrid
from mesa.time import RandomActivation

import random


def metastasis():
    value = random.randint(0, 50)
    if value == 0:
        return True


class TumorModel(Model):
    def __init__(self, width: int, height: int, num_agent: int):
        super().__init__()
        self.num_agent = num_agent
        self.schedule = RandomActivation(model=self)
        self.grid = MultiGrid(width=width, height=height,
                              torus=True)
        for agent_idx in range(num_agent):
            a = tumor.TumorCell(unique_id=agent_idx,
                                model=self, stem=True, metastatic=False)
            self.schedule.add(agent=a)
            x = int(height / 2)
            # x means the vertical axis here
            # (rotated cartesian because it's a matrix as if it was a picture)
            y = int(width / 2)
            self.grid.place_agent(agent=a,
                                  pos=(x, y))

    def step(self):
        self.schedule.step()
        for agent in self.schedule.agents:

            # metastasis
            metastasis_happening = metastasis()
            if agent.stem and metastasis_happening:
                a = tumor.TumorCell(unique_id=(self.num_agent + 1), model=self, stem=False, metastatic=True)
                self.num_agent += 1
                self.schedule.add(agent=a)
                cells_to_place = self.grid.get_neighborhood(
                    pos=agent.pos,
                    moore=False,
                    include_center=False,
                    radius=1)
                destination_cell = self.random.choice(cells_to_place)
                self.grid.place_agent(agent=a,
                                      pos=destination_cell)

            # division
            if agent.stem or agent.metastatic or (agent.age < 3 and agent.divide == 1):
                a = tumor.TumorCell(unique_id=(self.num_agent + 1),
                                    model=self, stem=False, metastatic=False)
                self.num_agent += 1
                self.schedule.add(agent=a)
                cells_to_place = self.grid.get_neighborhood(
                    pos=agent.pos,
                    moore=False,
                    include_center=False,
                    radius=1)
                destination_cell = self.random.choice(cells_to_place)
                self.grid.place_agent(agent=a,
                                      pos=destination_cell)

            # death
            if agent.age >= 7 and not (agent.stem or agent.metastatic):
                self.schedule.remove(agent=agent)
                self.grid.remove_agent(agent=agent)
