# Exceptions.py - All exceptions used by this library

class AlignmentException(Exception):
    """Base exception class for Alignments.py"""
    pass
    
class GeneNotFoundException(AlignmentException):
    """Raised when gene was not found in any of the alignments"""
    pass

class AlignmentParserException(Exception):
    """Base exception class for AlignmentParser.py"""
    pass

class FailedToParseException(AlignmentParserException):
    """Raised when failed to parse an alignment"""
    pass
    
class EmptyAlternativeTranscriptsException(AlignmentParserException):
    """
    Raised when calling Update() on the AlignmentRecord class when no
    alternative transcripts were found.
    """
    pass

class AlignmentRecordException(Exception):
    """Base exception class for RecordObjects.py"""
    pass
    
class AlternativeTranscriptRecordException(AlignmentRecordException):
    """Base exception class for alternative transcript records"""
    pass
   
class ZeroLengthedAlternativeTranscriptException(
        AlternativeTranscriptRecordException):
    """Raised when the spliced transcript has zero length"""
    pass