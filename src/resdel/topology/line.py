from typing import List, Optional


class Line:
    def __init__(self, name: str):
        self.raw: str
        self.tokens = self.raw.split()
        self.comment = None
        if ';' in self.tokens:
            self.comment = ' '.join(self.tokens[self.tokens.index(';') + 1:])
            self.tokens = self.tokens[:self.tokens.index(';')]

    def __repr__(self):
        return f"Line(raw={self.raw}, tokens={self.tokens}, comment={self.comment})"