import argparse
import nanome
from nanome_docking.Docking import Autodock4Docking, RhodiumDocking, SminaDocking


def main():
    parser = argparse.ArgumentParser(description='Parse Arguments to determine flavor of Docking to instantiate')
    parser.add_argument('algorithm', choices=['smina', 'autodock4'])

    args, _ = parser.parse_known_args()
    plugin_class = None
    algo = args.algorithm
    if algo == "smina":
        plugin_class = SminaDocking
    elif algo == "autodock4":
        plugin_class = Autodock4Docking
    elif algo == "rhodium":
        plugin_class = RhodiumDocking

    # Create the plugin, register Docking as the class to instantiate, and start listening
    plugin_name = f'Docking'
    description = f'Dock one or more ligands to a receptor using Smina or Autodock4.'
    category = "Docking"
    advanced_settings = True
    plugin = nanome.Plugin(plugin_name, description, category, advanced_settings)
    plugin.set_plugin_class(plugin_class)
    plugin.run()


if __name__ == "__main__":
    main()
