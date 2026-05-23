from pathlib import Path

class OutputPaths:
    def __init__(self, output_dir):
        self.root = Path(output_dir).resolve()

        self.root.mkdir(parents=True, exist_ok=True)

    def stage_dir(self, stage_name):
        path = self.root / stage_name
        path.mkdir(exist_ok=True)
        return path

    def topology_file(self, stage_name):
        return self.stage_dir(stage_name) / "system.top"
    
    def structure_PDBfile(self, stage_name):
        return self.stage_dir(stage_name) / "system.pdb"

    def structure_GROfile(self, stage_name):
        return self.stage_dir(stage_name) / "system.gro"
    
    def structure_POSREfile(self, stage_name):
        return self.stage_dir(stage_name) / "posre.itp"
