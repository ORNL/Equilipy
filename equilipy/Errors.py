class EquilibError(Exception):
    def __init__(self, message="Equilibrium calculation failed."):
        self.message = message
        super().__init__(self.message)
        
class DatabaseParsingError(Exception):
    def __init__(self, message="Failed to parse the database"):
        self.message = message
        super().__init__(self.message)
        
class InputConditionError(Exception):
    def __init__(self, message="Input condition is not recognized properly"):
        self.message = message
        super().__init__(self.message)
        
class PostProcessError(Exception):
    def __init__(self, message="Failed processing calcualtion result."):
        self.message = message
        super().__init__(self.message)
        
