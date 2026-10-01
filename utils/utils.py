import logging
import torch
import random
import numpy as np

def set_seed_everywhere(seed: int, device: str):
    """Set the seed for reproducibility."""
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if device == "cuda":
        torch.cuda.manual_seed_all(seed)

    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

def setup_logging(logging_path: str):
    """Sets up logging to a file and console."""
    root_logger = logging.getLogger("")
    root_logger.setLevel(logging.INFO)
    fmt = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

    # Add file handler (basicConfig is a no-op if handlers already exist)
    file_handler = logging.FileHandler(logging_path)
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(fmt)
    root_logger.addHandler(file_handler)

    # Add console handler if not already present
    if not any(isinstance(h, logging.StreamHandler) and not isinstance(h, logging.FileHandler) for h in root_logger.handlers):
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        console.setFormatter(fmt)
        root_logger.addHandler(console)

def to_tensor(array: np.array, device: str) -> torch.Tensor:
    """Convert np.array to torch.Tensor."""
    return torch.tensor(array, device=device)

def generate_causal_mask(size: int, device: str) -> torch.Tensor:
    """
    Generates the Causal Mask for input to Self-Attention. 
    Masked postions = float("-inf")
    Unmasked positions = float(0.0)

    :param size: Size of the (square) mask
    :return: A (size, size) mask
    """
    mask = (torch.triu(torch.ones(size, size, device=device)) == 1).transpose(0, 1)
    mask = mask.float().masked_fill(mask == 0, float("-inf")).masked_fill(mask == 1, float(0.0))
    return mask
