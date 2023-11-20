__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import pytest

import scine_molassembler as masm

def test_equal_molecules():
    graphA = "pGFhh6RhYQBhYwphcqNhbIKBCIELYmxygoEAgQFhc4KBCIELYXMBpGFhAGFjC2Fyo2Fsg4EJgQqBEGJscoOBAoEBgQBhc4OBEIEKgQlhcwKkYWEAYWMMYXKkYWyDgQeBD4EQY2xua4GiYXCCAQJjc2VxhQwPDg0QYmxyg4EAgQGBAmFzg4EHgQ+BEGFzA6RhYQBhYw1hcqRhbISBA4EEgQ6BEGNsbmuBomFwggIDY3NlcYUNDg8MEGJscoOCAAGBAoEDYXODggMEgQ6BEGFzBaRhYQBhYw5hcqRhbISBAYECgQ2BD2NsbmuBomFwggIDY3NlcYUODRAMD2JscoOCAAGBAoEDYXODggECgQ2BD2FzBaRhYQBhYw9hcqRhbISBBYEGgQyBDmNsbmuBomFwggIDY3NlcYUPDBANDmJscoOCAAGBA4ECYXODggUGgQ6BDGFzBaRhYQBhYxBhcqRhbISBAIELgQyBDWNsbmuBomFwggIDY3NlcYUQDA8ODWJscoSBAIEDgQGBAmFzhIEAgQ2BC4EMYXMFYWMPYWeiYUWRgwAQAIMBDgCDAg4AgwMNAIMEDQCDBQ8AgwYPAIMHDACDCAoAgwkLAIMKCwCDCxAAgwwPAIMMEACDDQ4Agw0QAIMODwBhWpEBAQEBAQEBAQEICAYHBgYGBmF2gwEBAA=="
    graphB = "pGFhg6RhYQBhYwdhcqNhbIOBBYEGgQlibHKDgQCBAoEBYXODgQWBCYEGYXMCpGFhAGFjCGFyo2FshIECgQOBBIEJYmxygoMAAQKBA2FzgoMCAwSBCWFzBaRhYQBhYwlhcqNhbISBAIEBgQeBCGJscoOCAAGBA4ECYXODggABgQiBB2FzBWFjD2FnomFFiYMACQCDAQkAgwIIAIMDCACDBAgAgwUHAIMGBwCDBwkAgwgJAGFaigEBAQEBAQgGBgZhdoMBAQA="
    graphC = "pGFhg6RhYQBhYwdhcqNhbIOBBYEGgQlibHKDgQCBAoEBYXODgQWBCYEGYXMCpGFhAGFjCGFyo2FshIECgQOBBIEJYmxygoMAAQKBA2FzgoMCAwSBCWFzBaRhYQBhYwlhcqNhbISBAIEBgQeBCGJscoOCAAGBA4ECYXODggABgQiBB2FzBWFjD2FnomFFiYMACQCDAQkAgwIIAIMDCACDBAgAgwUHAIMGBwCDBwkAgwgJAGFaigEBAQEBAQgGBgZhdoMBAgA=;pGFhh6RhYQBhYwphcqNhbIKBCIELYmxygoEAgQFhc4KBCIELYXMBpGFhAGFjC2Fyo2Fsg4EJgQqBEGJscoOBAoEBgQBhc4OBEIEKgQlhcwKkYWEAYWMMYXKkYWyDgQeBD4EQY2xua4GiYXCCAQJjc2VxhQwPDg0QYmxyg4EAgQGBAmFzg4EHgQ+BEGFzA6RhYQBhYw1hcqRhbISBA4EEgQ6BEGNsbmuBomFwggIDY3NlcYUNDg8MEGJscoOCAAGBAoEDYXODggMEgQ6BEGFzBaRhYQBhYw5hcqRhbISBAYECgQ2BD2NsbmuBomFwggIDY3NlcYUODRAMD2JscoOCAAGBAoEDYXODggECgQ2BD2FzBaRhYQBhYw9hcqRhbISBBYEGgQyBDmNsbmuBomFwggIDY3NlcYUPDBANDmJscoOCAAGBA4ECYXODggUGgQ6BDGFzBaRhYQBhYxBhcqRhbISBAIELgQyBDWNsbmuBomFwggIDY3NlcYUQDA8ODWJscoSBAIEDgQGBAmFzhIEAgQ2BC4EMYXMFYWMPYWeiYUWRgwAQAIMBDgCDAg4AgwMNAIMEDQCDBQ8AgwYPAIMHDACDCAoAgwkLAIMKCwCDCxAAgwwPAIMMEACDDQ4Agw0QAIMODwBhWpEBAQEBAQEBAQEICAYHBgYGBmF2gwECAA=="
    graphD = "pGFhg6RhYQBhYwdhcqNhbIOBBYEGgQlibHKDgQCBAoEBYXODgQWBCYEGYXMCpGFhAGFjCGFyo2FshIECgQOBBIEJYmxygoMAAQKBA2FzgoMCAwSBCWFzBaRhYQBhYwlhcqNhbISBAIEBgQeBCGJscoOCAAGBA4ECYXODggABgQiBB2FzBWFjD2FnomFFiYMACQCDAQkAgwIIAIMDCACDBAgAgwUHAIMGBwCDBwkAgwgJAGFaigEBAQEBAQgGBgZhdoMBAQA=;pGFhh6RhYQBhYwphcqNhbIKBCIELYmxygoEAgQFhc4KBCIELYXMBpGFhAGFjC2Fyo2Fsg4EJgQqBEGJscoOBAoEBgQBhc4OBEIEKgQlhcwKkYWEAYWMMYXKkYWyDgQeBD4EQY2xua4GiYXCCAQJjc2VxhQwPDg0QYmxyg4EAgQGBAmFzg4EHgQ+BEGFzA6RhYQBhYw1hcqRhbISBA4EEgQ6BEGNsbmuBomFwggIDY3NlcYUNDg8MEGJscoOCAAGBAoEDYXODggMEgQ6BEGFzBaRhYQBhYw5hcqRhbISBAYECgQ2BD2NsbmuBomFwggIDY3NlcYUODRAMD2JscoOCAAGBAoEDYXODggECgQ2BD2FzBaRhYQBhYw9hcqRhbISBBYEGgQyBDmNsbmuBomFwggIDY3NlcYUPDBANDmJscoOCAAGBA4ECYXODggUGgQ6BDGFzBaRhYQBhYxBhcqRhbISBAIELgQyBDWNsbmuBomFwggIDY3NlcYUQDA8ODWJscoSBAIEDgQGBAmFzhIEAgQ2BC4EMYXMFYWMPYWeiYUWRgwAQAIMBDgCDAg4AgwMNAIMEDQCDBQ8AgwYPAIMHDACDCAoAgwkLAIMKCwCDCxAAgwwPAIMMEACDDQ4Agw0QAIMODwBhWpEBAQEBAQEBAQEICAYHBgYGBmF2gwEBAA=="
    assert masm.JsonSerialization.equal_molecules(graphA, graphA)
    assert not masm.JsonSerialization.equal_molecules(graphA, graphB)
    assert not masm.JsonSerialization.equal_molecules(graphA, graphC)
    assert not masm.JsonSerialization.equal_molecules(graphA, graphD)
    assert not masm.JsonSerialization.equal_molecules(graphB, graphC)
    assert not masm.JsonSerialization.equal_molecules(graphB, graphD)
    assert masm.JsonSerialization.equal_molecules(graphC, graphC)
    assert masm.JsonSerialization.equal_molecules(graphC, graphD)
    assert masm.JsonSerialization.equal_molecules(graphD, graphC)

def test_equal_decision_lists():
    listA = "(52, 57, 63, 1):(54, 60, 65, 3)"
    listB = "(1, 7, 12, 1):(-165, -159, -154, 1)"
    listC = ";(70, 76, 81, 2)"
    listD = "(1, 7, 12, 1):(-165, -159, -154, 1);(70, 76, 81, 2)"
    listE = "(1, 7, 12, 1):(-165, -159, -154, 1);;(70, 76, 81, 2)"
    listF = "(173, 178, -177, 1)"
    listG = "(174, 179, -176, 1)"
    assert masm.JsonSerialization.equal_decision_lists(listA, listA)
    assert masm.JsonSerialization.equal_decision_lists(listB, listB)
    assert not masm.JsonSerialization.equal_decision_lists(listA, listB)
    assert not masm.JsonSerialization.equal_decision_lists(listB, listA)
    assert masm.JsonSerialization.equal_decision_lists(listC, listC)
    assert not masm.JsonSerialization.equal_decision_lists(listB, listD)
    assert masm.JsonSerialization.equal_decision_lists(listD, listD)
    assert not masm.JsonSerialization.equal_decision_lists(listE, listD)
    assert masm.JsonSerialization.equal_decision_lists(listF, listG)
