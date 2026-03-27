from typing import List, Optional


class Line:

    def __init__(self, line: str):

        self.raw = line.strip()
        self.tokens: Optional[List[str]] = None
        self.comment: Optional[str] = None

        self.is_section = False
        self.section_name: Optional[str] = None

        stripped = self.raw.strip()

        # empty line
        if not stripped:
            return

        # full comment
        if stripped.startswith(";"):
            return

        # section header
        if stripped.startswith("[") and stripped.endswith("]"):
            self.is_section = True
            self.section_name = stripped.strip("[] ").lower()
            return

        # inline comment
        if ";" in self.raw:
            data, comment = self.raw.split(";", 1)
            self.comment = comment.strip()
        else:
            data = self.raw

        tokens = data.split()

        if tokens:
            self.tokens = tokens

    def __repr__(self):

        return (
            f"Line(raw={self.raw}, tokens={self.tokens}, "
            f"comment={self.comment}, is_section={self.is_section}, "
            f"section_name={self.section_name})"
        )