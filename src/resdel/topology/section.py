from typing import List
from .line import Line


class Section:
    def __init__(self, name: str):
        self.name = name
        self.lines: List[Line] = []

    def add_line(self, line: Line):
        self.lines.append(line)

    def __repr__(self):
        return f"Section(name={self.name}, lines={self.lines})"