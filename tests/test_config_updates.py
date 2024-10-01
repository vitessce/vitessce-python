import pytest

from vitessce import (
    VitessceConfig,
    ViewType as cm,
)


@pytest.fixture
def vitessce_config():
    vc = VitessceConfig(schema_version="1.0.15")
    my_dataset = vc.add_dataset(name='My Dataset')
    vc.add_view(cm.SPATIAL, dataset=my_dataset)
    vc.add_view(cm.SCATTERPLOT, dataset=my_dataset)
    vc.add_view(cm.SCATTERPLOT, dataset=my_dataset)
    return vc


def test_get_views(vitessce_config):
    views = vitessce_config.get_views()
    assert len(views) == 3
    assert views[0].view["component"].lower() == "spatial"
    assert views[1].view["component"].lower() == "scatterplot"


def test_get_view_by_index(vitessce_config):
    view = vitessce_config.get_view_by_index(0)
    vc_dict = view.to_dict()
    assert vc_dict["component"] == "spatial"

    view = vitessce_config.get_view_by_index(1)
    vc_dict = view.to_dict()
    assert vc_dict["component"] == "scatterplot"

    with pytest.raises(IndexError):
        vitessce_config.get_view_by_index(5)


def test_get_view_by_component(vitessce_config):
    view = vitessce_config.get_first_view_by_type("spatial")
    vc_dict = view.to_dict()
    assert vc_dict["component"] == "spatial"

    view = vitessce_config.get_first_view_by_type("SCATTERPLOT")
    vc_dict = view.to_dict()
    assert vc_dict["component"] == "scatterplot"

    with pytest.raises(ValueError):
        vitessce_config.get_first_view_by_type("TEST")


def test_get_view_by_invalid_type(vitessce_config):
    with pytest.raises(TypeError):
        vitessce_config.get_first_view_by_type(3.5)


def test_remove_view_by_index(vitessce_config):
    removed_view = vitessce_config.remove_view_by_index(0)
    rv_dict = removed_view.to_dict()
    assert rv_dict["component"] == "spatial"
    assert len(vitessce_config.get_views()) == 2

    removed_view = vitessce_config.remove_view_by_index(1)
    rv_dict = removed_view.to_dict()
    assert rv_dict["component"] == "scatterplot"
    assert len(vitessce_config.get_views()) == 1

    with pytest.raises(IndexError):
        vitessce_config.remove_view_by_index(5)


def test_remove_view_by_component(vitessce_config):
    removed_view = vitessce_config.remove_first_view_by_type("spatial")
    rv_dict = removed_view.to_dict()
    assert rv_dict["component"] == "spatial"
    assert len(vitessce_config.get_views()) == 2

    with pytest.raises(ValueError):
        vitessce_config.remove_first_view_by_type("spatial")

    with pytest.raises(ValueError):
        vitessce_config.remove_first_view_by_type("TEST")


def test_remove_view_by_invalid_index(vitessce_config):
    with pytest.raises(TypeError):
        vitessce_config.remove_view_by_index(3.5)
