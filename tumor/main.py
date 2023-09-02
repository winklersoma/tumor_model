import matplotlib.pyplot as plt

import tumor


def main():
    width = 40
    height = 40
    model = tumor.TumorModel(width=width, height=height)
    time = 500  # number of times to run the model
    for t in range(time):
        model.step()
    model_data = model.datacollector.get_model_vars_dataframe()
    model_data.plot()
    plt.show()


if __name__ == "__main__":
    main()
