#!/usr/bin/env python
# coding: utf-8
#
# Reads the header from a GE a-Si Angio Detector
# Using version 8001 of the header from file:
#     c:\adept\core\DefaultImageInfoConfig.csv
#
#  Antonino Miceli
#  Thu Jan  4 13:46:31 CST 2007
#

# modifications by Jon Wright for style, pychecker and fabio
#
# Get ready for python3:
from __future__ import with_statement, print_function, division

__authors__ = ["Antonino Miceli" , "Jon Wright", "Jérôme Kieffer", "Harshad Paranjape"]
__date__ = "07/17/2015"
__status__ = "production"
__copyright__ = "2007 APS; 2010-2015 ESRF"
__licence__ = "GPL"


import numpy
import struct
import logging
import os
logger = logging.getLogger("GEimage")
from .fabioimage import fabioimage
from .fabioutils import next_filename, previous_filename

GE_HEADER_INFO = [
    # Name, length in bytes, format for struct (None means string)
    ('ImageFormat', 10, None),
    ('VersionOfStandardHeader', 2, '<H'),
    ('StandardHeaderSizeInBytes', 4, '<L'),
    ('VersionOfUserHeader', 2, '<H'),
    ('UserHeaderSizeInBytes', 4, '<L'),
    ('NumberOfFrames', 2, '<H'),
    ('NumberOfRowsInFrame', 2, '<H'),
    ('NumberOfColsInFrame', 2, '<H'),
    ('ImageDepthInBits', 2, '<H'),
    ('AcquisitionDate', 20, None),
    ('AcquisitionTime', 20, None),
    ('DUTID', 20, None),
    ('Operator', 50, None),
    ('DetectorSignature', 20, None),
    ('TestSystemName', 20, None),
    ('TestStationRevision', 20, None),
    ('CoreBundleRevision', 20, None),
    ('AcquisitionName', 40, None),
    ('AcquisitionParameterRevision', 20, None),
    ('OriginalNumberOfRows', 2, '<H'),
    ('OriginalNumberOfColumns', 2, '<H'),
    ('RowNumberUpperLeftPointArchiveROI', 2, '<H'),
    ('ColNumberUpperLeftPointArchiveROI', 2, '<H'),
    ('Swapped', 2, '<H'),
    ('Reordered', 2, '<H'),
    ('HorizontalFlipped', 2, '<H'),
    ('VerticalFlipped', 2, '<H'),
    ('WindowValueDesired', 2, '<H'),
    ('LevelValueDesired', 2, '<H'),
    ('AcquisitionMode', 2, '<H'),
    ('AcquisitionType', 2, '<H'),
    ('UserAcquisitionCoffFileName1', 100, None),
    ('UserAcquisitionCoffFileName2', 100, None),
    ('FramesBeforeExpose', 2, '<H'),
    ('FramesDuringExpose', 2, '<H'),
    ('FramesAfterExpose', 2, '<H'),
    ('IntervalBetweenFrames', 2, '<H'),
    ('ExposeTimeDelayInMicrosecs', 8, '<d'),
    ('TimeBetweenFramesInMicrosecs', 8, '<d'),
    ('FramesToSkipExpose', 2, '<H'),
    ('ExposureMode', 2, '<H'),
    ('PrepPresetTimeInMicrosecs', 8, '<d'),
    ('ExposePresetTimeInMicrosecs', 8, '<d'),
    ('AcquisitionFrameRateInFps', 4, '<f'),
    ('FOVSelect', 2, '<H'),
    ('ExpertMode', 2, '<H'),
    ('SetVCommon1', 8, '<d'),
    ('SetVCommon2', 8, '<d'),
    ('SetAREF', 8, '<d'),
    ('SetAREFTrim', 4, '<L'),
    ('SetSpareVoltageSource', 8, '<d'),
    ('SetCompensationVoltageSource', 8, '<d'),
    ('SetRowOffVoltage', 8, '<d'),
    ('SetRowOnVoltage', 8, '<d'),
    ('StoreCompensationVoltage', 4, '<L'),
    ('RampSelection', 2, '<H'),
    ('TimingMode', 2, '<H'),
    ('Bandwidth', 2, '<H'),
    ('ARCIntegrator', 2, '<H'),
    ('ARCPostIntegrator', 2, '<H'),
    ('NumberOfRows', 4, '<L'),
    ('RowEnable', 2, '<H'),
    ('EnableStretch', 2, '<H'),
    ('CompEnable', 2, '<H'),
    ('CompStretch', 2, '<H'),
    ('LeftEvenTristate', 2, '<H'),
    ('RightOddTristate', 2, '<H'),
    ('TestModeSelect', 4, '<L'),
    ('AnalogTestSource', 4, '<L'),
    ('VCommonSelect', 4, '<L'),
    ('DRCColumnSum', 4, '<L'),
    ('TestPatternFrameDelta', 4, '<L'),
    ('TestPatternRowDelta', 4, '<L'),
    ('TestPatternColumnDelta', 4, '<L'),
    ('DetectorHorizontalFlip', 2, '<H'),
    ('DetectorVerticalFlip', 2, '<H'),
    ('DFNAutoScrubOnOff', 2, '<H'),
    ('FiberChannelTimeOutInMicrosecs', 4, '<L'),
    ('DFNAutoScrubDelayInMicrosecs', 4, '<L'),
    ('StoreAECROI', 2, '<H'),
    ('TestPatternSaturationValue', 2, '<H'),
    ('TestPatternSeed', 4, '<L'),
    ('ExposureTimeInMillisecs', 4, '<f'),
    ('FrameRateInFps', 4, '<f'),
    ('kVp', 4, '<f'),
    ('mA', 4, '<f'),
    ('mAs', 4, '<f'),
    ('FocalSpotInMM', 4, '<f'),
    ('GeneratorType', 20, None),
    ('StrobeIntensityInFtL', 4, '<f'),
    ('NDFilterSelection', 2, '<H'),
    ('RefRegTemp1', 8, '<d'),
    ('RefRegTemp2', 8, '<d'),
    ('RefRegTemp3', 8, '<d'),
    ('Humidity1', 4, '<f'),
    ('Humidity2', 4, '<f'),
    ('DetectorControlTemp', 8, '<d'),
    ('DoseValueInmR', 8, '<d'),
    ('TargetLevelROIRow0', 2, '<H'),
    ('TargetLevelROICol0', 2, '<H'),
    ('TargetLevelROIRow1', 2, '<H'),
    ('TargetLevelROICol1', 2, '<H'),
    ('FrameNumberForTargetLevelROI', 2, '<H'),
    ('PercentRangeForTargetLevel', 2, '<H'),
    ('TargetValue', 2, '<H'),
    ('ComputedMedianValue', 2, '<H'),
    ('LoadZero', 2, '<H'),
    ('MaxLUTOut', 2, '<H'),
    ('MinLUTOut', 2, '<H'),
    ('MaxLinear', 2, '<H'),
    ('Reserved', 2, '<H'),
    ('ElectronsPerCount', 2, '<H'),
    ('ModeGain', 2, '<H'),
    ('TemperatureInDegC', 8, '<d'),
    ('LineRepaired', 2, '<H'),
    ('LineRepairFileName', 100, None),
    ('CurrentLongitudinalInMM', 4, '<f'),
    ('CurrentTransverseInMM', 4, '<f'),
    ('CurrentCircularInMM', 4, '<f'),
    ('CurrentFilterSelection', 4, '<L'),
    ('DisableScrubAck', 2, '<H'),
    ('ScanModeSelect', 2, '<H'),
    ('DetectorAppSwVersion', 20, None),
    ('DetectorNIOSVersion', 20, None),
    ('DetectorPeripheralSetVersion', 20, None),
    ('DetectorPhysicalAddress', 20, None),
    ('PowerDown', 2, '<H'),
    ('InitialVoltageLevel_VCOMMON', 8, '<d'),
    ('FinalVoltageLevel_VCOMMON', 8, '<d'),
    ('DmrCollimatorSpotSize', 10, None),
    ('DmrTrack', 5, None),
    ('DmrFilter', 5, None),
    ('FilterCarousel', 2, '<H'),
    ('Phantom', 20, None),
    ('SetEnableHighTime', 2, '<H'),
    ('SetEnableLowTime', 2, '<H'),
    ('SetCompHighTime', 2, '<H'),
    ('SetCompLowTime', 2, '<H'),
    ('SetSyncLowTime', 2, '<H'),
    ('SetConvertLowTime', 2, '<H'),
    ('SetSyncHighTime', 2, '<H'),
    ('SetEOLTime', 2, '<H'),
    ('SetRampOffsetTime', 2, '<H'),
    ('FOVStartingValue', 2, '<H'),
    ('ColumnBinning', 2, '<H'),
    ('RowBinning', 2, '<H'),
    ('BorderColumns64', 2, '<H'),
    ('BorderRows64', 2, '<H'),
    ('FETOffRows64', 2, '<H'),
    ('FOVStartColumn128', 2, '<H'),
    ('FOVStartRow128', 2, '<H'),
    ('NumberOfColumns128', 2, '<H'),
    ('NumberOfRows128', 2, '<H'),
    ('VFPAquisition', 2000, None),
    ('Comment', 200, None)
    ]

class GEimage(fabioimage):

    def __init__(self, filename=None,
                 NumberOfFrames=None, NumberOfRowsInFrame=None, NumberOfColsInFrame=None,
                 StandardHeaderSizeInBytes=None, UserHeaderSizeInBytes=None,
                 ImageDepthInBits=None):
        """
        Initialize GEimage with bare minimums. We don't want to rely on header
        for everything. We will accept params from user for number of frames, 
        number of rows/columns etc.
        """

        self._need_a_seek_to_read = True
        self.header = {}

        if filename is not None:
            self.filename = filename


        if StandardHeaderSizeInBytes is not None:
            self.header['StandardHeaderSizeInBytes'] = StandardHeaderSizeInBytes
        else:
            self.header['StandardHeaderSizeInBytes'] = 8192

        if UserHeaderSizeInBytes is not None:
            self.header['UserHeaderSizeInBytes'] = UserHeaderSizeInBytes
        else:
            self.header['UserHeaderSizeInBytes'] = 0

        if NumberOfRowsInFrame is not None:
            self.header['NumberOfRowsInFrame'] = NumberOfRowsInFrame
        else:
            self.header['NumberOfRowsInFrame'] = 2048

        if NumberOfColsInFrame is not None:
            self.header['NumberOfColsInFrame'] = NumberOfColsInFrame
        else:
            self.header['NumberOfColsInFrame'] = 2048

        if ImageDepthInBits is not None:
            self.header['ImageDepthInBits'] = ImageDepthInBits
        else:
            self.header['ImageDepthInBits'] = 16

        self.NumberOfBytesInFrame = self.header['NumberOfRowsInFrame'] * \
            self.header['NumberOfColsInFrame'] * \
            self.header['ImageDepthInBits'] // 8

        if NumberOfFrames is not None:
            self.header['NumberOfFrames'] = NumberOfFrames
            self.nframes = NumberOfFrames
        else:
            if self.filename is not None:
                self.header['NumberOfFrames'] = self.getNFrames()
                self.nframes = self.getNFrames()
            else:
                self.header['NumberOfFrames'] = 1
                self.nframes = 1


    def getNFrames(self):
        """
        Get number of frames in the file from file size
        """
        fileBytes = os.stat(self.filename).st_size
        nbytesFrame = self.NumberOfBytesInFrame
        nbytesHeader = self.header['StandardHeaderSizeInBytes']

        assert (fileBytes - nbytesHeader) % nbytesFrame == 0,\
            'file size not correct'
        nFrames = int((fileBytes - nbytesHeader) / nbytesFrame)
        if nFrames*nbytesFrame + nbytesHeader != fileBytes:
            raise RuntimeError, 'file size not correctly calculated'

        #nFrames = self.getNFramesFromBytes(fileBytes, self.header['StandardHeaderSizeInBytes'], self.NumberOfBytesInFrame)

        return nFrames

    def getNFramesFromBytes(fileBytes, nbytesHeader, nbytesFrame):
        """
        Calculate number of frames from total bytes, 
        header bytes and bytes in a frame
        """
        assert (fileBytes - nbytesHeader) % nbytesFrame == 0,\
            'file size not correct'
        nFrames = int((fileBytes - nbytesHeader) / nbytesFrame)
        if nFrames*nbytesFrame + nbytesHeader != fileBytes:
            raise RuntimeError, 'file size not correctly calculated'
        return nFrames


    def _readheader(self, infile):
        """ Read a GE image header """

        infile.seek(0)

        self.header = {}
        for name, nbytes, format in GE_HEADER_INFO:
            if format is None:
                self.header[ name ] = infile.read(nbytes)
            else:
                self.header[ name ] = struct.unpack(format,
                                                     infile.read(nbytes))[0]

    def read(self, fname=None, frame=None, readHeaderData=False):
        """
        Read in header into self.header and
        the data   into self.data
        """

        if fname is None and self.filename is None:
            raise RuntimeError, "No file to read"
        elif fname is None:
            fname = self.filename
        
        if readHeaderData:
            self.header = {}
            self.resetvals()
            self._readheader(infile)
            self.nframes = self.header['NumberOfFrames']

        infile = self._open(fname, "rb")
        self.sequencefilename = fname

        if frame is None:
            print('Reading %d frames from %s' % (self.header['NumberOfFrames'], fname))
            self._readstack(infile)
        else:
            print('Reading frame %d from %s' % (frame, fname))
            self._readframe(infile, frame)

        infile.close()

        return self

    def _makeframename(self):
        """ 
        Generate a frame name string by concatenating file name with
        fixed-width frame number
        """
        self.frameName = "%s$%04d" % (self.sequencefilename,
                                     self.currentFrameNum)

    def _readframe(self, filepointer, img_num):
        """
        Read a single frame from the image. The first image in the sequence is 0.
        This raises an exception if you give an invalid image; 
        otherwise fills in self.data
        """
        # Throw an exception for unusual frame number
        if(img_num > self.nframes):
            raise Exception("Image number out of bounds")
        elif(img_num < 0):
            raise Exception("Negative image number is not supported")
        # Determine at which point the frame data starts
        imgstart = self.header['StandardHeaderSizeInBytes'] + \
                   self.header['UserHeaderSizeInBytes'] + \
                   img_num * self.NumberOfBytesInFrame
        # Seek to the appropriate position from the beginning
        filepointer.seek(imgstart, 0)
        # Read the data in        
        data = numpy.fromfile(filepointer, count=self.header['NumberOfRowsInFrame']*self.header['NumberOfColsInFrame'], dtype=numpy.uint16).reshape(self.header['NumberOfRowsInFrame'], self.header['NumberOfColsInFrame'])
        # Take care on endianness
        if not numpy.little_endian:
            data.byteswap(True)

        data.shape = (self.header['NumberOfRowsInFrame'], self.header['NumberOfColsInFrame'])
        self.data = data
        self.dim2 , self.dim1 = self.data.shape
        self.currentFrameNum = int(img_num)
        self._makeframename()

    def _readstack(self, filepointer):
        """
        Read a stack of frames into a 3D data array
        """
        # Determine at which point the frame data starts
        imgstart = self.header['StandardHeaderSizeInBytes'] + \
                   self.header['UserHeaderSizeInBytes']

        # Seek to the appropriate position from the beginning
        filepointer.seek(imgstart, 0)
        # Read the data in        
        data = numpy.fromfile(filepointer, count=self.header['NumberOfRowsInFrame']*self.header['NumberOfColsInFrame']*self.header['NumberOfFrames'], dtype=numpy.uint16).reshape(self.header['NumberOfFrames'], self.header['NumberOfRowsInFrame'], self.header['NumberOfColsInFrame'])
        # Take care on endianness
        if not numpy.little_endian:
            data.byteswap(True)

        data.shape = (self.header['NumberOfFrames'], self.header['NumberOfRowsInFrame'], self.header['NumberOfColsInFrame'])
        self.data = data
        self.dim1 , self.dim2, self.dim3 = self.data.shape
        self.currentFrameNum = 0

    def write(self, fname=None, force_type=numpy.uint16):
        """ Not yet implemented"""
#        raise Exception("Write is not implemented")
        if fname is None and self.filename is None:
            raise RuntimeError, "There is no file to write to."
        elif fname is None:
            fname = self.filename

#        print("Writing data of shape [%d]" % (self.data.shape))
        self.data.tofile(fname)

    def getframe(self, img_num):
        """
        Returns a frame as a new fabioimage object
        """

        # Throw an exception for unusual frame number
        if(img_num > self.nframes):
            raise Exception("Image number out of bounds")
        elif(img_num < 0):
            raise Exception("Negative image number is not supported")

        # Do a deep copy of the header to make a new one
        newheader = {}
        for k in self.header.keys():
            newheader[k] = self.header[k]
        frame = GEimage(header=newheader)
        frame.nframes = self.nframes
        frame.sequencefilename = self.sequencefilename
        infile = frame._open(self.sequencefilename, "rb")
        frame._readframe(infile, img_num)
        infile.close()
        return frame

    def next(self):
        """
        Get the next image in a series as a fabio image
        """
        if self.currentFrameNum < (self.nframes - 1) and self.nframes > 1:
            return self.getframe(self.currentFrameNum + 1)
        else:
            newobj = GEimage()
            newobj.read(next_filename(
                self.sequencefilename))
            return newobj

    def previous(self):
        """
        Get the previous image in a series as a fabio image
        """
        if self.currentFrameNum > 0:
            return self.getframe(self.currentFrameNum - 1)
        else:
            newobj = GEimage()
            newobj.read(previous_filename(
                self.sequencefilename))
            return newobj


def demo():
    import sys, time

    if len(sys.argv) < 2:
        print("USAGE: GE_script.py <GEaSi_raw_image_file>")
        sys.exit()

    image_file = sys.argv[1]

    print("init read_GEaSi_data class and load header..")
    sequence1 = GEimage()
    sequence1.read(image_file)

    print("TimeBetweenFramesInMicrosecs = ")
    print(sequence1.header['TimeBetweenFramesInMicrosecs'])
    print("AcquisitionTime = ")
    print(sequence1.header['AcquisitionTime'])


    print("Mean = ", sequence1.data.ravel().mean())

    while 1:
        start = time.time()
        try:
            sequence1 = sequence1.next()
            print(sequence1.currentFrameNum, sequence1.data.ravel().mean(), \
                  time.time() - start)
        except Exception as  ex:
            raise ex




if __name__ == '__main__':
    demo()
