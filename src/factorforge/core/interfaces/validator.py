"""Validator interface"""

from abc import ABC, abstractmethod


class Validator(ABC):
    @abstractmethod
    def validate(self, data: str) -> bool:
        pass
