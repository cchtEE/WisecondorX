class ChromosomeMap:

    def __init__(self):
        self.max = 31
        self.memory = {i: str(i) for i in range(1, self.max + 1, 1)}
        self.memory_int = {i: str(i) for i in range(1, self.max + 1, 1)}

        self.memory[self.max - 1] = "X"
        self.memory[self.max] = "Y"

    def get_as_numeric_str(self, input_: int) -> str:
        return str(self.memory_int[input_])

    def get_as_numeric_int(self, input_: int) -> int:
        return int(self.memory_int[input_])

    def get_as_chr_name_str(self, input_: int) -> str:
        return str(self.memory[input_])

    def get_max(self):
        return self.max

    def get_x_index(self) -> int:
        return self.memory_int[self.max - 1]

    def get_y_index(self) -> int:
        return self.memory_int[self.max]
