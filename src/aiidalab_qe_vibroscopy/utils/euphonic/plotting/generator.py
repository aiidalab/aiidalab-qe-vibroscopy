from jinja2 import Environment, FileSystemLoader
import datetime
import os


def generate_from_template(model_state):
    # Load the Jinja environment and template
    env = Environment(loader=FileSystemLoader(os.path.dirname(__file__)))
    template = env.get_template("plot_heatmap.py.j2")

    # Define the context for the template
    context = {
        "generation_date": datetime.datetime.now().strftime("%Y-%m-%d"),
        "model_state": model_state,
    }

    # Render the template with the context
    output = template.render(context)

    return output


if __name__ == "__main__":
    # just for testing, we provide and empty metadata dictionary
    settings = {}
    script = generate_from_template(model_state=settings)
    print(script)
