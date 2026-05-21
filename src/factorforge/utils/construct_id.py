"""construct_id 생성 유틸리티."""
from datetime import datetime


def generate_construct_id() -> str:
    """CF-YYYYMMDD-HHMMSS 형식의 고유 construct ID 생성."""
    now = datetime.now()
    return f"CF-{now.strftime('%Y%m%d-%H%M%S')}"
