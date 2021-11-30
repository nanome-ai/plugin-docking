import asyncio
import os
import unittest
from unittest.mock import MagicMock, patch
from nanome.api.structure import Complex

from nanome_docking.Docking import Autodock4Docking, SminaDocking

fixtures_dir = os.path.join(os.path.dirname(__file__), 'fixtures')


class SminaDockingTestCase(unittest.TestCase):

    def setUp(self):
        receptor_sdf = f'{fixtures_dir}/5ceo_receptor.sdf'
        ligand_sdf = f'{fixtures_dir}/5ceo_ligand.sdf'
        self.receptor = Complex.io.from_sdf(path=receptor_sdf)
        self.ligand = Complex.io.from_sdf(path=ligand_sdf)

        self.plugin_smina = SminaDocking()
        self.plugin_smina.start()
        self.plugin_smina._network = MagicMock()

        self.plugin_autodock4 = Autodock4Docking()
        self.plugin_autodock4.start()
        self.plugin_autodock4._network = MagicMock()

    @patch('nanome.api.plugin_instance.PluginInstance.request_complexes')
    def test_smina_docking_complex(self, request_complexes_mock):
        print("Running Smina test")
        self.run_docking_test(self.plugin_smina, request_complexes_mock)

    @patch('nanome.api.plugin_instance.PluginInstance.request_complexes')
    def test_autodock4_docking(self, request_complexes_mock):
        print("Running autodock4 test")
        self.run_docking_test(self.plugin_autodock4, request_complexes_mock)

    def run_docking_test(self, plugin_instance, request_complexes_mock):
        """Test run_docking function on provided plugin_instance."""
        # Set future result for request_complexes mock
        fut = asyncio.Future()
        fut.set_result([self.receptor, self.ligand, self.ligand])
        request_complexes_mock.return_value = fut

        loop = asyncio.get_event_loop()
        mode_count = 2
        params = {
            'modes': mode_count,
            'exhaustiveness': 1,
            'autobox': 5,
        }
        result = loop.run_until_complete(
            plugin_instance.run_docking(self.receptor, [self.ligand], self.ligand, params)
        )
        comp = result[0]
        self.assertEqual(next(comp.molecules).conformer_count, mode_count)
