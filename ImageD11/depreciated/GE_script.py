#!/usr/bin/env python

#
# Reads the header from a GE a-Si Angio Detector
# Using version 8001 of the header from file:
#     c:\adept\core\DefaultImageInfoConfig.csv
#
#  Antonino Miceli
#  Thu Jan  4 13:46:31 CST 2007
#


# Minor modifications by Jon Wright for pychecker and
# adding some docstrings


import os,sys,Numeric,struct,mmap
from math import *

class read_GEaSi_data:
    """
    Class to read in an APS image format
    """
    def __init__(self,image_filename):
        """
        Reads the header 
        """
        self.image_filename = image_filename
        try:
            self.file = open(self.image_filename,"rb")
            self.size = os.path.getsize(self.image_filename)
            # print self.file.fileno(), self.size
            self.fp = mmap.mmap(self.file.fileno(), self.size,access=mmap.ACCESS_READ)

        except IOError, (errno, strerror):
            print "I/O error(%s): %s" % (errno, strerror)
            print "Error in opening file: ", self.image_filename
            raise

        self.fp.seek(0)
        # print self.fp.tell()

        # ADEPT
        self.ImageFormat = self.fp.read(10)

        # USHORT --> "=H"
        # ULONG  --> "=L"
        #   = means byte order is native

        self.VersionOfStandardHeader = self.fp.read(2)
        self.VersionOfStandardHeader = struct.unpack("=H",self.VersionOfStandardHeader)[0]

        self.StandardHeaderSizeInBytes = self.fp.read(4)
        self.StandardHeaderSizeInBytes = struct.unpack("=L",self.StandardHeaderSizeInBytes)[0]

        self.VersionOfUserHeader = self.fp.read(2)
        self.VersionOfUserHeader = struct.unpack("=H",self.VersionOfUserHeader)[0]

        self.UserHeaderSizeInBytes = self.fp.read(4)
        self.UserHeaderSizeInBytes = struct.unpack("=L",self.UserHeaderSizeInBytes)[0]

        self.NumberOfFrames = self.fp.read(2)
        self.NumberOfFrames =  struct.unpack("=H",self.NumberOfFrames)[0]

        self.NumberOfRowsInFrame  = self.fp.read(2)
        self.NumberOfRowsInFrame  =  struct.unpack("=H",self.NumberOfRowsInFrame)[0]

        self.NumberOfColsInFrame  = self.fp.read(2)
        self.NumberOfColsInFrame  =  struct.unpack("=H",self.NumberOfColsInFrame)[0]

        self.ImageDepthInBits  = self.fp.read(2)
        self.ImageDepthInBits  = struct.unpack("=H",self.ImageDepthInBits)[0]

        self.AcquisitionDate = self.fp.read(20)

        self.AcquisitionTime = self.fp.read(20)

        self.DUTID =  self.fp.read(20)

        self.Operator =  self.fp.read(50)

        self.DetectorSignature = self.fp.read(20)

        self.TestSystemName = self.fp.read(20)

        self.TestStationRevision = self.fp.read(20)

        self.CoreBundleRevision = self.fp.read(20)

        self.AcquisitionName = self.fp.read(40)

        self.AcquisitionParameterRevision = self.fp.read(20)

        self.OriginalNumberOfRows = self.fp.read(2)
        self.OriginalNumberOfRows = struct.unpack("=H",self.OriginalNumberOfRows)[0]

        self.OriginalNumberOfColumns = self.fp.read(2)
        self.OriginalNumberOfColumns = struct.unpack("=H",self.OriginalNumberOfColumns)[0]

        self.RowNumberUpperLeftPointArchiveROI = self.fp.read(2)
        self.RowNumberUpperLeftPointArchiveROI = struct.unpack("=H",self.RowNumberUpperLeftPointArchiveROI)[0]

        self.ColNumberUpperLeftPointArchiveROI = self.fp.read(2)
        self.ColNumberUpperLeftPointArchiveROI = struct.unpack("=H",self.ColNumberUpperLeftPointArchiveROI)[0]

        self.Swapped = self.fp.read(2) 
        self.Swapped = struct.unpack("=H",self.Swapped)[0]

        self.Reordered = self.fp.read(2) 
        self.Reordered = struct.unpack("=H",self.Reordered)[0]

        self.HorizontalFlipped = self.fp.read(2) 
        self.HorizontalFlipped = struct.unpack("=H",self.HorizontalFlipped)[0]

        self.VerticalFlipped = self.fp.read(2) 
        self.VerticalFlipped = struct.unpack("=H",self.VerticalFlipped)[0]

        self.WindowValueDesired = self.fp.read(2) 
        self.WindowValueDesired = struct.unpack("=H",self.WindowValueDesired)[0]

        self.LevelValueDesired = self.fp.read(2) 
        self.LevelValueDesired = struct.unpack("=H",self.LevelValueDesired)[0]

        self.AcquisitionMode = self.fp.read(2) 
        self.AcquisitionMode = struct.unpack("=H",self.AcquisitionMode)[0]

        self.AcquisitionType = self.fp.read(2) 
        self.AcquisitionType = struct.unpack("=H",self.AcquisitionType)[0]

        self.UserAcquisitionCoffFileName1 = self.fp.read(100) 
        self.UserAcquisitionCoffFileName2 = self.fp.read(100) 

        self.FramesBeforeExpose = self.fp.read(2) 
        self.FramesBeforeExpose = struct.unpack("=H",self.FramesBeforeExpose)[0]

        self.FramesDuringExpose = self.fp.read(2)  
        self.FramesDuringExpose = struct.unpack("=H",self.FramesDuringExpose)[0]

        self.FramesAfterExpose = self.fp.read(2) 
        self.FramesAfterExpose = struct.unpack("=H",self.FramesAfterExpose)[0]

        self.IntervalBetweenFrames = self.fp.read(2) 
        self.IntervalBetweenFrames = struct.unpack("=H",self.IntervalBetweenFrames)[0]

        self.ExposeTimeDelayInMicrosecs = self.fp.read(8) 
        self.ExposeTimeDelayInMicrosecs = struct.unpack("=d",self.ExposeTimeDelayInMicrosecs)[0]

        self.TimeBetweenFramesInMicrosecs = self.fp.read(8) 
        self.TimeBetweenFramesInMicrosecs = struct.unpack("=d",self.TimeBetweenFramesInMicrosecs)[0]

        self.FramesToSkipExpose = self.fp.read(2) 
        self.FramesToSkipExpose = struct.unpack("=H",self.FramesToSkipExpose)[0]

        # Rad --> ExposureMode = 1
        self.ExposureMode = self.fp.read(2) 
        self.ExposureMode = struct.unpack("=H",self.ExposureMode)[0]

        self.PrepPresetTimeInMicrosecs = self.fp.read(8) 
        self.PrepPresetTimeInMicrosecs = struct.unpack("=d",self.PrepPresetTimeInMicrosecs)[0]

        self.ExposePresetTimeInMicrosecs = self.fp.read(8) 
        self.ExposePresetTimeInMicrosecs = struct.unpack("=d",self.ExposePresetTimeInMicrosecs)[0]

        self.AcquisitionFrameRateInFps = self.fp.read(4) 
        self.AcquisitionFrameRateInFps = struct.unpack("=f",self.AcquisitionFrameRateInFps)[0]

        self.FOVSelect = self.fp.read(2)
        self.FOVSelect = struct.unpack("=H",self.FOVSelect)[0]

        self.ExpertMode = self.fp.read(2)
        self.ExpertMode = struct.unpack("=H",self.ExpertMode)[0]

        self.SetVCommon1 = self.fp.read(8)
        self.SetVCommon1 = struct.unpack("=d",self.SetVCommon1)[0]

        self.SetVCommon2 = self.fp.read(8)
        self.SetVCommon2 = struct.unpack("=d",self.SetVCommon2)[0]

        self.SetAREF = self.fp.read(8)
        self.SetAREF = struct.unpack("=d",self.SetAREF)[0]

        self.SetAREFTrim = self.fp.read(4)
        self.SetAREFTrim = struct.unpack("=L",self.SetAREFTrim)[0]

        self.SetSpareVoltageSource = self.fp.read(8)
        self.SetSpareVoltageSource = struct.unpack("=d",self.SetSpareVoltageSource)[0]

        self.SetCompensationVoltageSource = self.fp.read(8)
        self.SetCompensationVoltageSource = struct.unpack("=d",self.SetCompensationVoltageSource)[0]

        self.SetRowOffVoltage = self.fp.read(8)
        self.SetRowOffVoltage = struct.unpack("=d",self.SetRowOffVoltage)[0]

        self.SetRowOnVoltage = self.fp.read(8)
        self.SetRowOnVoltage = struct.unpack("=d",self.SetRowOnVoltage)[0]

        self.StoreCompensationVoltage = self.fp.read(4)
        self.StoreCompensationVoltage = struct.unpack("=L",self.StoreCompensationVoltage)[0]

        self.RampSelection = self.fp.read(2)
        self.RampSelection = struct.unpack("=H",self.RampSelection)[0]

        self.TimingMode = self.fp.read(2)
        self.TimingMode = struct.unpack("=H",self.TimingMode)[0]

        self.Bandwidth = self.fp.read(2)
        self.Bandwidth = struct.unpack("=H",self.Bandwidth)[0]

        self.ARCIntegrator = self.fp.read(2)
        self.ARCIntegrator = struct.unpack("=H",self.ARCIntegrator)[0]

        self.ARCPostIntegrator = self.fp.read(2)
        self.ARCPostIntegrator = struct.unpack("=H",self.ARCPostIntegrator)[0]

        self.NumberOfRows = self.fp.read(4)
        self.NumberOfRows = struct.unpack("=L",self.NumberOfRows)[0]

        self.RowEnable = self.fp.read(2)
        self.RowEnable = struct.unpack("=H",self.RowEnable)[0]

        self.EnableStretch = self.fp.read(2)
        self.EnableStretch = struct.unpack("=H",self.EnableStretch)[0]

        self.CompEnable = self.fp.read(2)
        self.CompEnable = struct.unpack("=H",self.CompEnable)[0]

        self.CompStretch = self.fp.read(2)
        self.CompStretch = struct.unpack("=H",self.CompStretch)[0]

        self.LeftEvenTristate = self.fp.read(2)
        self.LeftEvenTristate = struct.unpack("=H",self.LeftEvenTristate)[0]

        self.RightOddTristate = self.fp.read(2)
        self.RightOddTristate = struct.unpack("=H",self.RightOddTristate)[0]

        self.TestModeSelect = self.fp.read(4)
        self.TestModeSelect = struct.unpack("=L",self.TestModeSelect)[0]

        self.AnalogTestSource = self.fp.read(4)
        self.AnalogTestSource = struct.unpack("=L",self.AnalogTestSource)[0]

        self.VCommonSelect = self.fp.read(4)
        self.VCommonSelect = struct.unpack("=L",self.VCommonSelect)[0]

        self.DRCColumnSum = self.fp.read(4)
        self.DRCColumnSum = struct.unpack("=L",self.DRCColumnSum)[0]

        self.TestPatternFrameDelta = self.fp.read(4)
        self.TestPatternFrameDelta = struct.unpack("=L",self.TestPatternFrameDelta)[0]

        self.TestPatternRowDelta = self.fp.read(4)
        self.TestPatternRowDelta = struct.unpack("=L",self.TestPatternRowDelta)[0]

        self.TestPatternColumnDelta = self.fp.read(4)
        self.TestPatternColumnDelta = struct.unpack("=L",self.TestPatternColumnDelta)[0]

        self.DetectorHorizontalFlip = self.fp.read(2)
        self.DetectorHorizontalFlip = struct.unpack("=H",self.DetectorHorizontalFlip)[0]

        self.DetectorVerticalFlip = self.fp.read(2)
        self.DetectorVerticalFlip = struct.unpack("=H",self.DetectorVerticalFlip)[0]

        self.DFNAutoScrubOnOff = self.fp.read(2)
        self.DFNAutoScrubOnOff = struct.unpack("=H",self.DFNAutoScrubOnOff)[0]

        self.FiberChannelTimeOutInMicrosecs = self.fp.read(4)
        self.FiberChannelTimeOutInMicrosecs = struct.unpack("=L",self.FiberChannelTimeOutInMicrosecs)[0]

        self.DFNAutoScrubDelayInMicrosecs = self.fp.read(4)
        self.DFNAutoScrubDelayInMicrosecs = struct.unpack("=L",self.DFNAutoScrubDelayInMicrosecs)[0]

        self.StoreAECROI = self.fp.read(2)
        self.StoreAECROI = struct.unpack("=H",self.StoreAECROI)[0]

        self.TestPatternSaturationValue = self.fp.read(2)
        self.TestPatternSaturationValue = struct.unpack("=H",self.TestPatternSaturationValue)[0]

        self.TestPatternSeed = self.fp.read(4)
        self.TestPatternSeed = struct.unpack("=L",self.TestPatternSeed)[0]

        self.ExposureTimeInMillisecs = self.fp.read(4) 
        self.ExposureTimeInMillisecs = struct.unpack("=f",self.ExposureTimeInMillisecs)[0]

        self.FrameRateInFps = self.fp.read(4) 
        self.FrameRateInFps = struct.unpack("=f",self.FrameRateInFps)[0]

        self.kVp = self.fp.read(4) 
        self.kVp = struct.unpack("=f",self.kVp)[0]

        self.mA = self.fp.read(4) 
        self.mA = struct.unpack("=f",self.mA)[0]

        self.mAs = self.fp.read(4) 
        self.mAs = struct.unpack("=f",self.mAs)[0]

        self.FocalSpotInMM = self.fp.read(4) 
        self.FocalSpotInMM = struct.unpack("=f",self.FocalSpotInMM)[0]

        self.GeneratorType = self.fp.read(20)

        self.StrobeIntensityInFtL = self.fp.read(4) 
        self.StrobeIntensityInFtL = struct.unpack("=f",self.StrobeIntensityInFtL)[0]

        self.NDFilterSelection = self.fp.read(2) 
        self.NDFilterSelection = struct.unpack("=H",self.NDFilterSelection)[0]

        self.RefRegTemp1 = self.fp.read(8) 
        self.RefRegTemp1 = struct.unpack("=d",self.RefRegTemp1)[0]

        self.RefRegTemp2 = self.fp.read(8) 
        self.RefRegTemp2 = struct.unpack("=d",self.RefRegTemp2)[0]

        self.RefRegTemp3 = self.fp.read(8) 
        self.RefRegTemp3 = struct.unpack("=d",self.RefRegTemp3)[0]

        self.Humidity1 = self.fp.read(4) 
        self.Humidity1 = struct.unpack("=f",self.Humidity1)[0]

        self.Humidity2 = self.fp.read(4) 
        self.Humidity2 = struct.unpack("=f",self.Humidity2)[0]

        self.DetectorControlTemp = self.fp.read(8) 
        self.DetectorControlTemp = struct.unpack("=d",self.DetectorControlTemp)[0]

        self.DoseValueInmR = self.fp.read(8) 
        self.DoseValueInmR = struct.unpack("=d",self.DoseValueInmR)[0]

        self.TargetLevelROIRow0 = self.fp.read(2)
        self.TargetLevelROIRow0 = struct.unpack("=H",self.TargetLevelROIRow0)[0]

        self.TargetLevelROICol0 = self.fp.read(2)
        self.TargetLevelROICol0 = struct.unpack("=H",self.TargetLevelROICol0)[0]

        self.TargetLevelROIRow1 = self.fp.read(2)
        self.TargetLevelROIRow1 = struct.unpack("=H",self.TargetLevelROIRow1)[0]

        self.TargetLevelROICol1 = self.fp.read(2)
        self.TargetLevelROICol1 = struct.unpack("=H",self.TargetLevelROICol1)[0]

        self.FrameNumberForTargetLevelROI = self.fp.read(2)
        self.FrameNumberForTargetLevelROI = struct.unpack("=H",self.FrameNumberForTargetLevelROI)[0]

        self.PercentRangeForTargetLevel = self.fp.read(2)
        self.PercentRangeForTargetLevel = struct.unpack("=H",self.PercentRangeForTargetLevel)[0]

        self.TargetValue = self.fp.read(2)
        self.TargetValue = struct.unpack("=H",self.TargetValue)[0]

        self.ComputedMedianValue = self.fp.read(2)
        self.ComputedMedianValue = struct.unpack("=H",self.ComputedMedianValue)[0]

        self.LoadZero = self.fp.read(2)
        self.LoadZero = struct.unpack("=H",self.LoadZero)[0]

        self.MaxLUTOut = self.fp.read(2)
        self.MaxLUTOut = struct.unpack("=H",self.MaxLUTOut)[0]

        self.MinLUTOut = self.fp.read(2)
        self.MinLUTOut = struct.unpack("=H",self.MinLUTOut)[0]

        self.MaxLinear = self.fp.read(2)
        self.MaxLinear = struct.unpack("=H",self.MaxLinear)[0]

        self.Reserved = self.fp.read(2)
        self.Reserved = struct.unpack("=H",self.Reserved)[0]

        self.ElectronsPerCount = self.fp.read(2)
        self.ElectronsPerCount = struct.unpack("=H",self.ElectronsPerCount)[0]

        self.ModeGain = self.fp.read(2)
        self.ModeGain = struct.unpack("=H",self.ModeGain)[0]

        self.TemperatureInDegC = self.fp.read(8)
        self.TemperatureInDegC = struct.unpack("=d",self.TemperatureInDegC)[0]

        self.LineRepaired = self.fp.read(2)
        self.LineRepaired = struct.unpack("=H",self.LineRepaired)[0]

        self.LineRepairFileName = self.fp.read(100)

        self.CurrentLongitudinalInMM = self.fp.read(4)
        self.CurrentLongitudinalInMM = struct.unpack("=f",self.CurrentLongitudinalInMM)[0]

        self.CurrentTransverseInMM = self.fp.read(4)
        self.CurrentTransverseInMM = struct.unpack("=f",self.CurrentTransverseInMM)[0]

        self.CurrentCircularInMM = self.fp.read(4)
        self.CurrentCircularInMM = struct.unpack("=f",self.CurrentCircularInMM)[0]

        self.CurrentFilterSelection = self.fp.read(4)
        self.CurrentFilterSelection = struct.unpack("=L",self.CurrentFilterSelection)[0]

        self.DisableScrubAck = self.fp.read(2)
        self.DisableScrubAck = struct.unpack("=H",self.DisableScrubAck)[0]

        self.ScanModeSelect = self.fp.read(2)
        self.ScanModeSelect = struct.unpack("=H",self.ScanModeSelect)[0]

        self.DetectorAppSwVersion = self.fp.read(20)	

        self.DetectorNIOSVersion = self.fp.read(20)	

        self.DetectorPeripheralSetVersion = self.fp.read(20)	

        self.DetectorPhysicalAddress	 = self.fp.read(20)

        self.PowerDown = self.fp.read(2)
        self.PowerDown = struct.unpack("=H",self.PowerDown)[0]

        self.InitialVoltageLevel_VCOMMON = self.fp.read(8)
        self.InitialVoltageLevel_VCOMMON = struct.unpack("=d",self.InitialVoltageLevel_VCOMMON)[0]

        self.FinalVoltageLevel_VCOMMON = self.fp.read(8)
        self.FinalVoltageLevel_VCOMMON = struct.unpack("=d",self.FinalVoltageLevel_VCOMMON)[0]

        self.DmrCollimatorSpotSize	 = self.fp.read(10)

        self.DmrTrack	 = self.fp.read(5)

        self.DmrFilter	 = self.fp.read(5)

        self.FilterCarousel = self.fp.read(2)
        self.FilterCarousel = struct.unpack("=H",self.FilterCarousel)[0]

        self.Phantom	 = self.fp.read(20)

        self.SetEnableHighTime = self.fp.read(2)
        self.SetEnableHighTime = struct.unpack("=H",self.SetEnableHighTime)[0]

        self.SetEnableLowTime = self.fp.read(2)
        self.SetEnableLowTime = struct.unpack("=H",self.SetEnableLowTime)[0]

        self.SetCompHighTime = self.fp.read(2)
        self.SetCompHighTime = struct.unpack("=H",self.SetCompHighTime)[0]

        self.SetCompLowTime = self.fp.read(2)
        self.SetCompLowTime = struct.unpack("=H",self.SetCompLowTime)[0]

        self.SetSyncLowTime = self.fp.read(2)
        self.SetSyncLowTime = struct.unpack("=H",self.SetSyncLowTime)[0]

        self.SetConvertLowTime = self.fp.read(2)
        self.SetConvertLowTime = struct.unpack("=H",self.SetConvertLowTime)[0]

        self.SetSyncHighTime = self.fp.read(2)
        self.SetSyncHighTime = struct.unpack("=H",self.SetSyncHighTime)[0]

        self.SetEOLTime = self.fp.read(2)
        self.SetEOLTime = struct.unpack("=H",self.SetEOLTime)[0]

        self.SetRampOffsetTime = self.fp.read(2)
        self.SetRampOffsetTime = struct.unpack("=H",self.SetRampOffsetTime)[0]

        self.FOVStartingValue = self.fp.read(2)
        self.FOVStartingValue = struct.unpack("=H",self.FOVStartingValue)[0]

        self.ColumnBinning = self.fp.read(2)
        self.ColumnBinning = struct.unpack("=H",self.ColumnBinning)[0]

        self.RowBinning = self.fp.read(2)
        self.RowBinning = struct.unpack("=H",self.RowBinning)[0]

        self.BorderColumns64 = self.fp.read(2)
        self.BorderColumns64 = struct.unpack("=H",self.BorderColumns64)[0]

        self.BorderRows64 = self.fp.read(2)
        self.BorderRows64 = struct.unpack("=H",self.BorderRows64)[0]

        self.FETOffRows64 = self.fp.read(2)
        self.FETOffRows64 = struct.unpack("=H",self.FETOffRows64)[0]

        self.FOVStartColumn128 = self.fp.read(2)
        self.FOVStartColumn128 = struct.unpack("=H",self.FOVStartColumn128)[0]

        self.FOVStartRow128 = self.fp.read(2)
        self.FOVStartRow128 = struct.unpack("=H",self.FOVStartRow128)[0]

        self.NumberOfColumns128 = self.fp.read(2)
        self.NumberOfColumns128 = struct.unpack("=H",self.NumberOfColumns128)[0]

        self.NumberOfRows128 = self.fp.read(2)
        self.NumberOfRows128 = struct.unpack("=H",self.NumberOfRows128)[0]

        self.VFPAquisition	 = self.fp.read(2000)

        self.Comment	 = self.fp.read(200)

        # print fp.tell()
        # fp.seek(8192)


    def load_image_from_seq(self,img_num):
        # Load only one image from the sequence
        #    Note: the first image in the sequence is 1  (not 0 !!)
        # Returns -1 if you give an invalid image #, otherwise returns a Numeric array

        img_num = int(img_num)

        if(img_num > self.NumberOfFrames or img_num < 1):
            return(-1)

        # Go to the beginning of the file
        self.fp.seek(0)
        self.fp.seek(self.StandardHeaderSizeInBytes + self.UserHeaderSizeInBytes +
                     (img_num-1)*self.NumberOfRowsInFrame * self.NumberOfColsInFrame * (self.ImageDepthInBits/8)) 


        image_tmp = self.fp.read(self.NumberOfRowsInFrame * self.NumberOfColsInFrame * (self.ImageDepthInBits/8) )
        # unpack binary into 16-bit unsigned ints, little-endian
        format = "<%dH" % (self.NumberOfRowsInFrame * self.NumberOfColsInFrame)
        image_tmp = struct.unpack(format,image_tmp)
        image_tmp = Numeric.array(image_tmp)
        return(image_tmp)


    def close(self):
        self.fp.close()


if __name__ == '__main__':

    if len(sys.argv)<2:
        print "USAGE: read_GEaSi_data.py <GEaSi_raw_image_file>"
        sys.exit()

    image_file = sys.argv[1]

    print "init read_GEaSi_data class and load header.."
    sequence1 = read_GEaSi_data(image_file)    

    print "TimeBetweenFramesInMicrosecs = ",
    print sequence1.TimeBetweenFramesInMicrosecs
    print "AcquisitionTime = ",
    print sequence1.AcquisitionTime

    print "load_image_data..."
    image3 = sequence1.load_image_from_seq(1)

    if(image3 == -1):
        print "You gave an invalid image number... too large or negative or zero....Try again"
        sys.exit(1)
    import scipy
    print "Median = ",
    print scipy.median(image3)
    print "Mean = ",
    print scipy.mean(image3)
    print "Std = ",
    print scipy.std(image3)
    print "\n"

    image3 = Numeric.reshape(image3, (sequence1.NumberOfRowsInFrame,sequence1.NumberOfColsInFrame)) 
    print len(image3), image3[0][0], image3[1][1],image3[2047][2047]

   
