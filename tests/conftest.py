import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--ebtel_idl_path", action="store", default=None, help="Path to EBTEL IDL code"
    )
    parser.addoption(
        "--show_plots", action="store_true", default=False,
    )


@pytest.fixture
def ebtel_idl_path(request):
    return request.config.getoption("--ebtel_idl_path")


@pytest.fixture
def show_plots(request):
    return request.config.getoption("--show_plots")
