class Example:
    def __init__(self, parameter):
        # EITHER
        # variant a to initialize var_1
        self.var_1 = self.initialize_var_1_variant_a(parameter)
        # OR
        # variant b to initialize var_1
        self.initialize_var_1_variant_b(parameter)
        # OR something else
        # ...

    def initialize_var_1_variant_a(self, parameter):
        # complex calculations, var_1 = f(parameter)
        result_of_complex_calculations = 123
        return result_of_complex_calculations

    def initialize_var_1_variant_b(self, parameter):
        # complex calculations, var_1 = f(parameter)
        result_of_complex_calculations = 123
        self.var_1 = result_of_complex_calculations


example_instance = Example("some_parameter")
print(example_instance.var_1)