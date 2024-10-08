<?xml version="1.0" encoding="ISO-8859-1"?>

<xsl:stylesheet version="1.0"
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

    <xsl:output method="xml" indent="yes"/>

    <xsl:variable name="phaseOversampling">
        <xsl:choose>
            <xsl:when test="not(siemens/IRIS/DERIVED/phaseOversampling)">0</xsl:when>
            <xsl:otherwise>
                <xsl:value-of select="siemens/IRIS/DERIVED/phaseOversampling"/>
            </xsl:otherwise>
        </xsl:choose>
    </xsl:variable>

    <xsl:variable name="sliceOversampling">
        <xsl:choose>
            <xsl:when test="not(siemens/MEAS/sKSpace/dSliceOversamplingForDialog)">0</xsl:when>
            <xsl:otherwise>
                <xsl:value-of select="siemens/MEAS/sKSpace/dSliceOversamplingForDialog"/>
            </xsl:otherwise>
        </xsl:choose>
    </xsl:variable>

    <xsl:variable name="partialFourierPhase">
        <xsl:choose>
<!-- Kelvin: These value are incorrect, at least in VD/VE
            <xsl:when test="siemens/MEAS/sKSpace/ucPhasePartialFourier = 1">0.5</xsl:when>
            <xsl:when test="siemens/MEAS/sKSpace/ucPhasePartialFourier = 2">0.75</xsl:when>
            <xsl:when test="siemens/MEAS/sKSpace/ucPhasePartialFourier = 4">0.875</xsl:when>
-->
            <xsl:when test="siemens/MEAS/sKSpace/ucPhasePartialFourier = 1">0.5</xsl:when>
            <xsl:when test="siemens/MEAS/sKSpace/ucPhasePartialFourier = 2">0.625</xsl:when>
            <xsl:when test="siemens/MEAS/sKSpace/ucPhasePartialFourier = 4">0.75</xsl:when>
            <xsl:when test="siemens/MEAS/sKSpace/ucPhasePartialFourier = 8">0.875</xsl:when>
            <xsl:otherwise>1.0</xsl:otherwise>
        </xsl:choose>
    </xsl:variable>

    <xsl:variable name="numberOfContrasts">
        <xsl:value-of select="siemens/MEAS/lContrasts"/>
    </xsl:variable>

    <xsl:variable name="studyID">
        <xsl:value-of select="substring(siemens/IRIS/RECOMPOSE/StudyLOID, 6)"/>
    </xsl:variable>

    <xsl:variable name="patientID">
        <xsl:value-of select="substring(siemens/IRIS/RECOMPOSE/PatientLOID, 6)"/>
    </xsl:variable>

    <xsl:variable name="strSeperator">_</xsl:variable>

    <xsl:template match="/">
        <ismrmrdHeader xsi:schemaLocation="http://www.ismrm.org/ISMRMRD ismrmrd.xsd"
                xmlns="http://www.ismrm.org/ISMRMRD"
                xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
                xmlns:xs="http://www.w3.org/2001/XMLSchema">

            <subjectInformation>
                <patientName>
                    <xsl:value-of select="$patientID"/>
                </patientName>
                <xsl:if test="siemens/YAPS/flUsedPatientWeight > 0">
                    <patientWeight_kg>
                        <xsl:value-of select="siemens/YAPS/flUsedPatientWeight"/>
                    </patientWeight_kg>
                </xsl:if>
                <patientID>
                    <xsl:value-of select="$patientID"/>
                </patientID>
                <patientGender>
                    <xsl:choose>
                        <xsl:when test="siemens/DICOM/lPatientSex = 1">F</xsl:when>
                        <xsl:when test="siemens/DICOM/lPatientSex = 2">M</xsl:when>
                        <xsl:otherwise>O</xsl:otherwise>
                    </xsl:choose>
                </patientGender>
            </subjectInformation>

            <studyInformation>
                <studyInstanceUID>
                    <xsl:value-of select="$studyID" />
                </studyInstanceUID>
            </studyInformation>

            <measurementInformation>
                <measurementID>
                    <xsl:value-of select="concat(string(siemens/DICOM/DeviceSerialNumber), $strSeperator, $patientID, $strSeperator, $studyID, $strSeperator, string(siemens/HEADER/MeasUID))"/>
                </measurementID>
                <patientPosition>
                    <xsl:value-of select="siemens/YAPS/tPatientPosition"/>
                </patientPosition>
                <protocolName>
                    <xsl:value-of select="siemens/MEAS/tProtocolName"/>
                </protocolName>

                <xsl:if test="siemens/YAPS/ReconMeasDependencies/RFMap > 0">
                    <measurementDependency>
                        <dependencyType>RFMap</dependencyType>
                        <measurementID>
                            <xsl:value-of select="concat(string(siemens/DICOM/DeviceSerialNumber), $strSeperator, $patientID, $strSeperator, $studyID, $strSeperator, string(siemens/YAPS/ReconMeasDependencies/RFMap))"/>
                        </measurementID>
                    </measurementDependency>
                </xsl:if>

                <xsl:if test="siemens/YAPS/ReconMeasDependencies/SenMap > 0">
                    <measurementDependency>
                        <dependencyType>SenMap</dependencyType>
                        <measurementID>
                            <xsl:value-of select="concat(string(siemens/DICOM/DeviceSerialNumber), $strSeperator, $patientID, $strSeperator, $studyID, $strSeperator, string(siemens/YAPS/ReconMeasDependencies/SenMap))"/>
                        </measurementID>
                    </measurementDependency>
                </xsl:if>

                <xsl:if test="siemens/YAPS/ReconMeasDependencies/Noise > 0">
                    <measurementDependency>
                        <dependencyType>Noise</dependencyType>
                        <measurementID>
                            <xsl:value-of select="concat(string(siemens/DICOM/DeviceSerialNumber), $strSeperator, $patientID, $strSeperator, $studyID, $strSeperator, string(siemens/YAPS/ReconMeasDependencies/Noise))"/>
                        </measurementID>
                    </measurementDependency>
                </xsl:if>

                <frameOfReferenceUID>
                    <xsl:value-of select="siemens/YAPS/tFrameOfReference" />
                </frameOfReferenceUID>

            </measurementInformation>

            <acquisitionSystemInformation>
                <systemVendor>
                    <xsl:value-of select="siemens/DICOM/Manufacturer"/>
                </systemVendor>
                <systemModel>
                    <xsl:value-of select="siemens/DICOM/ManufacturersModelName"/>
                </systemModel>
                <systemFieldStrength_T>
                    <xsl:value-of select="siemens/YAPS/flMagneticFieldStrength"/>
                </systemFieldStrength_T>
                <relativeReceiverNoiseBandwidth>0.793</relativeReceiverNoiseBandwidth>
                <receiverChannels>
                    <xsl:value-of select="siemens/YAPS/iMaxNoOfRxChannels" />
                </receiverChannels>

                <!-- Coil Labels -->
                <xsl:choose>
                    <!-- VD line with dual density -->
                    <xsl:when test="siemens/MEAS/asCoilSelectMeas/ADC/lADCChannelConnected">
                        <xsl:variable name="NumberOfSelectedCoils">
                            <xsl:value-of select="count(siemens/MEAS/asCoilSelectMeas/Select/lElementSelected[text() = '1'])" />
                        </xsl:variable>
                        <xsl:for-each select="siemens/MEAS/asCoilSelectMeas/ADC/lADCChannelConnected[position() >= 1  and not(position() > $NumberOfSelectedCoils)]">
                            <xsl:sort data-type="number"
                                      select="." />
                            <xsl:variable name="CurADC"
                                          select="."/>
                            <xsl:variable name="CurADCIndex"
                                          select="position()" />
                            <xsl:for-each select="../lADCChannelConnected[position() >= 1  and not(position() > $NumberOfSelectedCoils)]">
                                <xsl:if test="$CurADC = .">
                                    <xsl:variable name="CurCoil" select="position()"/>
                                    <xsl:variable name="CurCoilID" select="../../ID/tCoilID[$CurCoil]"/>
                                    <xsl:variable name="CurCoilElement" select="../../Elem/tElement[$CurCoil]"/>
                                    <xsl:variable name="CurCoilCopyID" select="../../Coil/lCoilCopy[$CurCoil]"/>
                                    <coilLabel>
                                        <coilNumber>
                                            <xsl:value-of select="number(../lADCChannelConnected[$CurADCIndex])"/>
                                        </coilNumber>
                                        <coilName>
                                            <xsl:value-of select="$CurCoilID"/>:<xsl:value-of select="string($CurCoilCopyID)"/>:<xsl:value-of select="$CurCoilElement"/>
                                        </coilName>
                                    </coilLabel>
                                </xsl:if>
                            </xsl:for-each>
                        </xsl:for-each>
                    </xsl:when>
                    <xsl:otherwise>
                        <!-- This is probably VB -->
                        <xsl:for-each select="siemens/MEAS/asCoilSelectMeas/ID/tCoilID">
                            <xsl:variable name="CurCoil"
                                          select="position()"/>
                            <coilLabel>
                                <coilNumber>
                                    <xsl:value-of select="$CurCoil -1"/>
                                </coilNumber>
                                <coilName>
                                    <xsl:value-of select="."/>:<xsl:value-of select="../../Elem/tElement[$CurCoil]"/>
                                </coilName>
                            </coilLabel>
                        </xsl:for-each>
                    </xsl:otherwise>
                </xsl:choose>

                <institutionName>
                    <xsl:value-of select="siemens/DICOM/InstitutionName" />
                </institutionName>
            </acquisitionSystemInformation>

            <experimentalConditions>
                <H1resonanceFrequency_Hz>
                    <xsl:value-of select="siemens/DICOM/lFrequency"/>
                </H1resonanceFrequency_Hz>
            </experimentalConditions>
            <encoding>
                <trajectory>
                    <xsl:choose>
                        <xsl:when test="siemens/MEAS/sKSpace/ucTrajectory = 1">cartesian</xsl:when>
                        <xsl:when test="siemens/MEAS/sKSpace/ucTrajectory = 2">radial</xsl:when>
                        <xsl:when test="siemens/MEAS/sKSpace/ucTrajectory = 4">spiral</xsl:when>
                        <xsl:when test="siemens/MEAS/sKSpace/ucTrajectory = 8">propellor</xsl:when>
                        <xsl:otherwise>other</xsl:otherwise>
                    </xsl:choose>
                </trajectory>

                <xsl:if test="siemens/MEAS/sKSpace/ucTrajectory = 4">
                    <trajectoryDescription>
                        <identifier>HargreavesVDS2000</identifier>
                        <userParameterLong>
                            <name>interleaves</name>
                            <value>
                                <xsl:value-of select="siemens/MEAS/sKSpace/lRadialViews" />
                            </value>
                        </userParameterLong>
                        <userParameterLong>
                            <name>fov_coefficients</name>
                            <value>1</value>
                        </userParameterLong>
                        <userParameterLong>
                            <name>SamplingTime_ns</name>
                            <value>
                                <xsl:value-of select="siemens/MEAS/sWipMemBlock/alFree[57]" />
                            </value>
                        </userParameterLong>
                        <userParameterDouble>
                            <name>MaxGradient_G_per_cm</name>
                            <value>
                                <xsl:value-of select="siemens/MEAS/sWipMemBlock/adFree[7]" />
                            </value>
                        </userParameterDouble>
                        <userParameterDouble>
                            <name>MaxSlewRate_G_per_cm_per_s</name>
                            <value>
                                <xsl:value-of select="siemens/MEAS/sWipMemBlock/adFree[8]" />
                            </value>
                        </userParameterDouble>
                        <userParameterDouble>
                            <name>FOVCoeff_1_cm</name>
                            <value>
                                <xsl:value-of select="siemens/MEAS/sWipMemBlock/adFree[10]" />
                            </value>
                        </userParameterDouble>
                        <userParameterDouble>
                            <name>krmax_per_cm</name>
                            <value>
                                <xsl:value-of select="siemens/MEAS/sWipMemBlock/adFree[9]" />
                            </value>
                        </userParameterDouble>
                        <comment>Using spiral design by Brian Hargreaves (http://mrsrl.stanford.edu/~brian/vdspiral/)</comment>
                    </trajectoryDescription>
                </xsl:if>

                <xsl:if test="siemens/YAPS/alRegridRampupTime > 0">
                    <xsl:if test="siemens/YAPS/alRegridRampdownTime > 0">
                        <trajectoryDescription>
                            <identifier>ConventionalEPI</identifier>
                            <userParameterLong>
                                <name>etl</name>
                                <value>
                                    <xsl:value-of select="siemens/MEAS/sFastImaging/lEPIFactor"/>
                                </value>
                            </userParameterLong>
                            <userParameterLong>
                                <name>numberOfNavigators</name>
                                <value>3</value>
                            </userParameterLong>
                            <userParameterLong>
                                <name>rampUpTime</name>
                                <value>
                                    <xsl:value-of select="siemens/YAPS/alRegridRampupTime"/>
                                </value>
                            </userParameterLong>
                            <userParameterLong>
                                <name>rampDownTime</name>
                                <value>
                                    <xsl:value-of select="siemens/YAPS/alRegridRampdownTime"/>
                                </value>
                            </userParameterLong>
                            <userParameterLong>
                                <name>flatTopTime</name>
                                <value>
                                    <xsl:value-of select="siemens/YAPS/alRegridFlattopTime"/>
                                </value>
                            </userParameterLong>
                            <userParameterLong>
                                <name>echoSpacing</name>
                                <value>
                                    <xsl:value-of select="siemens/YAPS/lEchoSpacing"/>
                                </value>
                            </userParameterLong>
                            <userParameterLong>
                                <name>acqDelayTime</name>
                                <value>
                                    <xsl:value-of select="siemens/YAPS/alRegridDelaySamplesTime"/>
                                </value>
                            </userParameterLong>
                            <userParameterLong>
                                <name>numSamples</name>
                                <value>
                                    <xsl:value-of select="siemens/YAPS/alRegridDestSamples"/>
                                </value>
                            </userParameterLong>
                            <userParameterDouble>
                                <name>dwellTime</name>
                                <value>
                                    <xsl:value-of select="siemens/MEAS/sRXSPEC/alDwellTime div 1000.0"/>
                                </value>
                            </userParameterDouble>
                            <comment>Conventional 2D EPI sequence</comment>
                        </trajectoryDescription>
                    </xsl:if>
                </xsl:if>

                <encodedSpace>
                    <matrixSize>

                        <xsl:choose>
                            <xsl:when test="siemens/MEAS/sKSpace/ucTrajectory = 1">
                                <x>
                                    <xsl:value-of select="siemens/YAPS/iNoOfFourierColumns"/>
                                </x>
                            </xsl:when>
                            <xsl:otherwise>
                                <x>
                                    <xsl:value-of select="siemens/IRIS/DERIVED/imageColumns"/>
                                </x>
                            </xsl:otherwise>
                        </xsl:choose>

                        <xsl:choose>
                            <xsl:when test="siemens/MEAS/sKSpace/uc2DInterpolation" >
                                <xsl:choose>
                                    <xsl:when test="siemens/MEAS/sKSpace/uc2DInterpolation = 1">
                                        <y>
                                            <xsl:value-of select="floor(siemens/YAPS/iPEFTLength div 2)"/>
                                        </y>
                                    </xsl:when>
                                    <xsl:otherwise>
                                        <y>
                                            <xsl:value-of select="siemens/YAPS/iPEFTLength"/>
                                        </y>
                                    </xsl:otherwise>
                                </xsl:choose>
                            </xsl:when>
                            <xsl:otherwise>
                                <y>
                                    <xsl:value-of select="siemens/YAPS/iPEFTLength"/>
                                </y>
                            </xsl:otherwise>
                        </xsl:choose>

                        <xsl:choose>
                            <xsl:when test="not(siemens/YAPS/iNoOfFourierPartitions) or (siemens/YAPS/i3DFTLength = 1)">
                                <z>1</z>
                            </xsl:when>
                            <xsl:otherwise>
                                <z>
                                    <xsl:value-of select="siemens/YAPS/i3DFTLength"/>
                                </z>
                            </xsl:otherwise>
                        </xsl:choose>
                    </matrixSize>

                    <fieldOfView_mm>
                        <xsl:choose>
                            <xsl:when test="siemens/MEAS/sKSpace/ucTrajectory = 1">
                                <x>
                                    <xsl:value-of select="siemens/MEAS/sSliceArray/asSlice/s0/dReadoutFOV * siemens/YAPS/flReadoutOSFactor"/>
                                </x>
                            </xsl:when>
                            <xsl:otherwise>
                                <x>
                                    <xsl:value-of select="siemens/MEAS/sSliceArray/asSlice/s0/dReadoutFOV"/>
                                </x>
                            </xsl:otherwise>
                        </xsl:choose>
                        <y>
                            <xsl:value-of select="siemens/MEAS/sSliceArray/asSlice/s0/dPhaseFOV * (1+$phaseOversampling)"/>
                        </y>
                        <z>
                            <xsl:value-of select="siemens/MEAS/sSliceArray/asSlice/s0/dThickness * (1+$sliceOversampling)"/>
                        </z>
                    </fieldOfView_mm>
                </encodedSpace>
                <reconSpace>
                    <matrixSize>
                        <x>
                            <xsl:value-of select="siemens/IRIS/DERIVED/imageColumns"/>
                        </x>
                        <y>
                            <xsl:value-of select="siemens/IRIS/DERIVED/imageLines"/>
                        </y>
                        <xsl:choose>
                            <xsl:when test="siemens/YAPS/i3DFTLength = 1">
                                <z>1</z>
                            </xsl:when>
                            <xsl:otherwise>
                                <z>
                                    <xsl:value-of select="siemens/MEAS/sKSpace/lImagesPerSlab"/>
                                </z>
                            </xsl:otherwise>
                        </xsl:choose>
                    </matrixSize>
                    <fieldOfView_mm>
                        <x>
                            <xsl:value-of select="siemens/MEAS/sSliceArray/asSlice/s0/dReadoutFOV"/>
                        </x>
                        <y>
                            <xsl:value-of select="siemens/MEAS/sSliceArray/asSlice/s0/dPhaseFOV"/>
                        </y>
                        <z>
                            <xsl:value-of select="siemens/MEAS/sSliceArray/asSlice/s0/dThickness"/>
                        </z>
                    </fieldOfView_mm>
                </reconSpace>
                <encodingLimits>
                    <kspace_encoding_step_1>
                        <minimum>0</minimum>
                        <maximum>
                            <xsl:value-of select="siemens/YAPS/iNoOfFourierLines - 1"/>
                        </maximum>
                        <center>
<!-- Kelvin: Different logic because the shorter half of partial Fourier is acquired first -->
<!--                             <xsl:value-of select="floor(siemens/MEAS/sKSpace/lPhaseEncodingLines div 2)"/> -->
<!--                             <xsl:value-of select="ceiling((siemens/YAPS/iNoOfFourierLines - 1) div 2 * $partialFourierPhase)"/> -->
                             <xsl:value-of select="ceiling((siemens/YAPS/iNoOfFourierLines - 1) * (1 - (0.5 div $partialFourierPhase)) )"/>
                        </center>
                    </kspace_encoding_step_1>
                    <kspace_encoding_step_2>
                        <minimum>0</minimum>
                        <xsl:choose>
                            <xsl:when test="not(siemens/YAPS/iNoOfFourierPartitions) or (siemens/YAPS/i3DFTLength = 1)">
                                <maximum>0</maximum>
                                <center>0</center>
                            </xsl:when>
                            <xsl:otherwise>
                                <maximum>
                                    <xsl:value-of select="siemens/YAPS/iNoOfFourierPartitions - 1"/>
                                </maximum>
                                <xsl:choose>
                                    <xsl:when test="siemens/MEAS/sKSpace/ucTrajectory = 1">
                                        <xsl:choose>
                                            <xsl:when test="siemens/MEAS/sPat/lAccelFact3D">
                                                <xsl:choose>
                                                    <xsl:when test="not(siemens/MEAS/sPat/lAccelFact3D) > 1">
                                                        <center>
                                                            <xsl:value-of select="floor(siemens/MEAS/sKSpace/lPartitions div 2) - (siemens/YAPS/lPartitions - siemens/YAPS/iNoOfFourierPartitions)"/>
                                                        </center>
                                                    </xsl:when>
                                                    <xsl:otherwise>
                                                        <xsl:choose>
                                                            <xsl:when test="(siemens/MEAS/sKSpace/lPartitions - siemens/YAPS/iNoOfFourierPartitions) > siemens/MEAS/sPat/lAccelFact3D">
                                                                <center>
                                                                    <xsl:value-of select="floor(siemens/MEAS/sKSpace/lPartitions div 2) - (siemens/MEAS/sKSpace/lPartitions - siemens/YAPS/iNoOfFourierPartitions)"/>
                                                                </center>
                                                            </xsl:when>
                                                            <xsl:otherwise>
                                                                <center>
                                                                    <xsl:value-of select="floor(siemens/MEAS/sKSpace/lPartitions div 2)"/>
                                                                </center>
                                                            </xsl:otherwise>
                                                        </xsl:choose>
                                                    </xsl:otherwise>
                                                </xsl:choose>
                                            </xsl:when>
                                            <xsl:otherwise>
                                                <center>
                                                    <xsl:value-of select="floor(siemens/MEAS/sKSpace/lPartitions div 2) - (siemens/MEAS/sKSpace/lPartitions - siemens/YAPS/iNoOfFourierPartitions)"/>
                                                </center>
                                            </xsl:otherwise>
                                        </xsl:choose>
                                    </xsl:when>
                                    <xsl:otherwise>
                                        <center>0</center>
                                    </xsl:otherwise>
                                </xsl:choose>
                            </xsl:otherwise>
                        </xsl:choose>
                    </kspace_encoding_step_2>
                    <slice>
                        <minimum>0</minimum>
                        <maximum>
                            <xsl:value-of select="siemens/MEAS/sSliceArray/lSize - 1"/>
                        </maximum>
                        <center>0</center>
                    </slice>
                    <set>
                        <minimum>0</minimum>
                        <maximum>
                            <xsl:choose>
                                <xsl:when test="siemens/YAPS/iNSet">
                                    <xsl:value-of select="siemens/YAPS/iNSet - 1"/>
                                </xsl:when>
                                <xsl:otherwise>0</xsl:otherwise>
                            </xsl:choose>
                        </maximum>
                        <center>0</center>
                    </set>
                    <phase>
                        <minimum>0</minimum>
                        <maximum>
                            <xsl:choose>
                                <xsl:when test="siemens/MEAS/sPhysioImaging/lPhases">
                                    <xsl:value-of select="siemens/MEAS/sPhysioImaging/lPhases - 1"/>
                                </xsl:when>
                                <xsl:otherwise>0</xsl:otherwise>
                            </xsl:choose>
                        </maximum>
                        <center>0</center>
                    </phase>
                    <repetition>
                        <minimum>0</minimum>
                        <maximum>
                            <xsl:choose>
                                <xsl:when test="siemens/MEAS/lRepetitions">
                                    <xsl:value-of select="siemens/MEAS/lRepetitions"/>
                                </xsl:when>
                                <xsl:otherwise>0</xsl:otherwise>
                            </xsl:choose>
                        </maximum>
                        <center>0</center>
                    </repetition>
                    <segment>
                        <minimum>0</minimum>
                        <maximum>
                            <xsl:choose>
                                <xsl:when test="siemens/MEAS/sFastImaging/ucSegmentationMode" >
                                    <xsl:choose>
                                        <xsl:when test="siemens/MEAS/sFastImaging/ucSegmentationMode = 2">
                                            <xsl:choose>
                                                <xsl:when test="siemens/MEAS/sFastImaging/lShots">
                                                    <xsl:value-of select="siemens/MEAS/sFastImaging/lShots - 1"/>
                                                </xsl:when>
                                                <xsl:otherwise>0</xsl:otherwise>
                                            </xsl:choose>
                                        </xsl:when>
                                        <xsl:when test="siemens/MEAS/sFastImaging/ucSegmentationMode = 1">
                                            <xsl:choose>
                                                <xsl:when test="siemens/MEAS/sFastImaging/lSegments &gt; 1">
                                                    <xsl:value-of select="ceiling((siemens/YAPS/iNoOfFourierPartitions * siemens/YAPS/iNoOfFourierLines) div siemens/MEAS/sFastImaging/lSegments)"/>
                                                </xsl:when>
                                                <xsl:otherwise>0</xsl:otherwise>
                                            </xsl:choose>
                                        </xsl:when>
                                        <xsl:otherwise>0</xsl:otherwise>
                                    </xsl:choose>
                                </xsl:when>
                                <xsl:otherwise>0</xsl:otherwise>
                            </xsl:choose>
                        </maximum>
                        <center>0</center>
                    </segment>
                    <contrast>
                        <minimum>0</minimum>
                        <maximum>
                            <xsl:choose>
                                <xsl:when test="siemens/MEAS/lContrasts">
                                    <xsl:value-of select="siemens/MEAS/lContrasts - 1"/>
                                </xsl:when>
                                <xsl:otherwise>0</xsl:otherwise>
                            </xsl:choose>
                        </maximum>
                        <center>0</center>
                    </contrast>
                    <average>
                        <minimum>0</minimum>
                        <maximum>
                            <xsl:choose>
                                <xsl:when test="siemens/MEAS/lAverages">
                                    <xsl:value-of select="siemens/MEAS/lAverages - 1"/>
                                </xsl:when>
                                <xsl:otherwise>0</xsl:otherwise>
                            </xsl:choose>
                        </maximum>
                        <center>0</center>
                    </average>
                </encodingLimits>
        <parallelImaging>
          <accelerationFactor>
                    <kspace_encoding_step_1>
              <xsl:choose>
            <xsl:when test="not(siemens/MEAS/sPat/lAccelFactPE)">1</xsl:when>
            <xsl:otherwise>
              <xsl:value-of select="(siemens/MEAS/sPat/lAccelFactPE)"/>
            </xsl:otherwise>
              </xsl:choose>
                    </kspace_encoding_step_1>
                    <kspace_encoding_step_2>
              <xsl:choose>
            <xsl:when test="not(siemens/MEAS/sPat/lAccelFact3D)">1</xsl:when>
            <xsl:otherwise>
              <xsl:value-of select="(siemens/MEAS/sPat/lAccelFact3D)"/>
            </xsl:otherwise>
              </xsl:choose>
                    </kspace_encoding_step_2>
          </accelerationFactor>
          <calibrationMode>
                    <xsl:choose>
              <xsl:when test="siemens/MEAS/sPat/ucRefScanMode = 1">other</xsl:when>
              <xsl:when test="siemens/MEAS/sPat/ucRefScanMode = 2">embedded</xsl:when>
              <xsl:when test="siemens/MEAS/sPat/ucRefScanMode = 4">separate</xsl:when>
              <xsl:when test="siemens/MEAS/sPat/ucRefScanMode = 8">separate</xsl:when>
              <xsl:when test="siemens/MEAS/sPat/ucRefScanMode = 16">interleaved</xsl:when>
              <xsl:when test="siemens/MEAS/sPat/ucRefScanMode = 32">interleaved</xsl:when>
              <xsl:when test="siemens/MEAS/sPat/ucRefScanMode = 64">interleaved</xsl:when>
              <xsl:otherwise>other</xsl:otherwise>
                    </xsl:choose>
          </calibrationMode>
          <xsl:if test="(siemens/MEAS/sPat/ucRefScanMode = 1) or (siemens/MEAS/sPat/ucRefScanMode = 16) or (siemens/MEAS/sPat/ucRefScanMode = 32) or (siemens/MEAS/sPat/ucRefScanMode = 64)">
                    <interleavingDimension>
              <xsl:choose>
            <xsl:when test="siemens/MEAS/sPat/ucRefScanMode = 16">average</xsl:when>
            <xsl:when test="siemens/MEAS/sPat/ucRefScanMode = 32">repetition</xsl:when>
            <xsl:when test="siemens/MEAS/sPat/ucRefScanMode = 64">phase</xsl:when>
            <xsl:otherwise>other</xsl:otherwise>
              </xsl:choose>
                    </interleavingDimension>
          </xsl:if>
        </parallelImaging>
            </encoding>

            <!-- Encoding space 1: High-contrast SASHA lines acquired with acceleration rate 3 -->
            <encoding>
                <trajectory>
                    <xsl:choose>
                        <xsl:when test="siemens/MEAS/sKSpace/ucTrajectory = 1">cartesian</xsl:when>
                        <xsl:when test="siemens/MEAS/sKSpace/ucTrajectory = 2">radial</xsl:when>
                        <xsl:when test="siemens/MEAS/sKSpace/ucTrajectory = 4">spiral</xsl:when>
                        <xsl:when test="siemens/MEAS/sKSpace/ucTrajectory = 8">propellor</xsl:when>
                        <xsl:otherwise>other</xsl:otherwise>
                    </xsl:choose>
                </trajectory>

                <xsl:if test="siemens/MEAS/sKSpace/ucTrajectory = 4">
                    <trajectoryDescription>
                        <identifier>HargreavesVDS2000</identifier>
                        <userParameterLong>
                            <name>interleaves</name>
                            <value>
                                <xsl:value-of select="siemens/MEAS/sKSpace/lRadialViews" />
                            </value>
                        </userParameterLong>
                        <userParameterLong>
                            <name>fov_coefficients</name>
                            <value>1</value>
                        </userParameterLong>
                        <userParameterLong>
                            <name>SamplingTime_ns</name>
                            <value>
                                <xsl:value-of select="siemens/MEAS/sWipMemBlock/alFree[57]" />
                            </value>
                        </userParameterLong>
                        <userParameterDouble>
                            <name>MaxGradient_G_per_cm</name>
                            <value>
                                <xsl:value-of select="siemens/MEAS/sWipMemBlock/adFree[7]" />
                            </value>
                        </userParameterDouble>
                        <userParameterDouble>
                            <name>MaxSlewRate_G_per_cm_per_s</name>
                            <value>
                                <xsl:value-of select="siemens/MEAS/sWipMemBlock/adFree[8]" />
                            </value>
                        </userParameterDouble>
                        <userParameterDouble>
                            <name>FOVCoeff_1_cm</name>
                            <value>
                                <xsl:value-of select="siemens/MEAS/sWipMemBlock/adFree[10]" />
                            </value>
                        </userParameterDouble>
                        <userParameterDouble>
                            <name>krmax_per_cm</name>
                            <value>
                                <xsl:value-of select="siemens/MEAS/sWipMemBlock/adFree[9]" />
                            </value>
                        </userParameterDouble>
                        <comment>Using spiral design by Brian Hargreaves (http://mrsrl.stanford.edu/~brian/vdspiral/)</comment>
                    </trajectoryDescription>
                </xsl:if>

                <xsl:if test="siemens/YAPS/alRegridRampupTime > 0">
                    <xsl:if test="siemens/YAPS/alRegridRampdownTime > 0">
                        <trajectoryDescription>
                            <identifier>ConventionalEPI</identifier>
                            <userParameterLong>
                                <name>etl</name>
                                <value>
                                    <xsl:value-of select="siemens/MEAS/sFastImaging/lEPIFactor"/>
                                </value>
                            </userParameterLong>
                            <userParameterLong>
                                <name>numberOfNavigators</name>
                                <value>3</value>
                            </userParameterLong>
                            <userParameterLong>
                                <name>rampUpTime</name>
                                <value>
                                    <xsl:value-of select="siemens/YAPS/alRegridRampupTime"/>
                                </value>
                            </userParameterLong>
                            <userParameterLong>
                                <name>rampDownTime</name>
                                <value>
                                    <xsl:value-of select="siemens/YAPS/alRegridRampdownTime"/>
                                </value>
                            </userParameterLong>
                            <userParameterLong>
                                <name>flatTopTime</name>
                                <value>
                                    <xsl:value-of select="siemens/YAPS/alRegridFlattopTime"/>
                                </value>
                            </userParameterLong>
                            <userParameterLong>
                                <name>echoSpacing</name>
                                <value>
                                    <xsl:value-of select="siemens/YAPS/lEchoSpacing"/>
                                </value>
                            </userParameterLong>
                            <userParameterLong>
                                <name>acqDelayTime</name>
                                <value>
                                    <xsl:value-of select="siemens/YAPS/alRegridDelaySamplesTime"/>
                                </value>
                            </userParameterLong>
                            <userParameterLong>
                                <name>numSamples</name>
                                <value>
                                    <xsl:value-of select="siemens/YAPS/alRegridDestSamples"/>
                                </value>
                            </userParameterLong>
                            <userParameterDouble>
                                <name>dwellTime</name>
                                <value>
                                    <xsl:value-of select="siemens/MEAS/sRXSPEC/alDwellTime div 1000.0"/>
                                </value>
                            </userParameterDouble>
                            <comment>Conventional 2D EPI sequence</comment>
                        </trajectoryDescription>
                    </xsl:if>
                </xsl:if>

                <encodedSpace>
                    <matrixSize>

                        <xsl:choose>
                            <xsl:when test="siemens/MEAS/sKSpace/ucTrajectory = 1">
                                <x>
                                    <xsl:value-of select="siemens/YAPS/iNoOfFourierColumns"/>
                                </x>
                            </xsl:when>
                            <xsl:otherwise>
                                <x>
                                    <xsl:value-of select="siemens/IRIS/DERIVED/imageColumns"/>
                                </x>
                            </xsl:otherwise>
                        </xsl:choose>

                        <xsl:choose>
                            <xsl:when test="siemens/MEAS/sKSpace/uc2DInterpolation" >
                                <xsl:choose>
                                    <xsl:when test="siemens/MEAS/sKSpace/uc2DInterpolation = 1">
                                        <y>
                                            <xsl:value-of select="floor(siemens/YAPS/iPEFTLength div 2)"/>
                                        </y>
                                    </xsl:when>
                                    <xsl:otherwise>
                                        <y>
                                            <xsl:value-of select="siemens/YAPS/iPEFTLength"/>
                                        </y>
                                    </xsl:otherwise>
                                </xsl:choose>
                            </xsl:when>
                            <xsl:otherwise>
                                <y>
                                    <xsl:value-of select="siemens/YAPS/iPEFTLength"/>
                                </y>
                            </xsl:otherwise>
                        </xsl:choose>

                        <xsl:choose>
                            <xsl:when test="not(siemens/YAPS/iNoOfFourierPartitions) or (siemens/YAPS/i3DFTLength = 1)">
                                <z>1</z>
                            </xsl:when>
                            <xsl:otherwise>
                                <z>
                                    <xsl:value-of select="siemens/YAPS/i3DFTLength"/>
                                </z>
                            </xsl:otherwise>
                        </xsl:choose>
                    </matrixSize>

                    <fieldOfView_mm>
                        <xsl:choose>
                            <xsl:when test="siemens/MEAS/sKSpace/ucTrajectory = 1">
                                <x>
                                    <xsl:value-of select="siemens/MEAS/sSliceArray/asSlice/s0/dReadoutFOV * siemens/YAPS/flReadoutOSFactor"/>
                                </x>
                            </xsl:when>
                            <xsl:otherwise>
                                <x>
                                    <xsl:value-of select="siemens/MEAS/sSliceArray/asSlice/s0/dReadoutFOV"/>
                                </x>
                            </xsl:otherwise>
                        </xsl:choose>
                        <y>
                            <xsl:value-of select="siemens/MEAS/sSliceArray/asSlice/s0/dPhaseFOV * (1+$phaseOversampling)"/>
                        </y>
                        <z>
                            <xsl:value-of select="siemens/MEAS/sSliceArray/asSlice/s0/dThickness * (1+$sliceOversampling)"/>
                        </z>
                    </fieldOfView_mm>
                </encodedSpace>
                <reconSpace>
                    <matrixSize>
                        <x>
                            <xsl:value-of select="siemens/IRIS/DERIVED/imageColumns"/>
                        </x>
                        <y>
                            <xsl:value-of select="siemens/IRIS/DERIVED/imageLines"/>
                        </y>
                        <xsl:choose>
                            <xsl:when test="siemens/YAPS/i3DFTLength = 1">
                                <z>1</z>
                            </xsl:when>
                            <xsl:otherwise>
                                <z>
                                    <xsl:value-of select="siemens/MEAS/sKSpace/lImagesPerSlab"/>
                                </z>
                            </xsl:otherwise>
                        </xsl:choose>
                    </matrixSize>
                    <fieldOfView_mm>
                        <x>
                            <xsl:value-of select="siemens/MEAS/sSliceArray/asSlice/s0/dReadoutFOV"/>
                        </x>
                        <y>
                            <xsl:value-of select="siemens/MEAS/sSliceArray/asSlice/s0/dPhaseFOV"/>
                        </y>
                        <z>
                            <xsl:value-of select="siemens/MEAS/sSliceArray/asSlice/s0/dThickness"/>
                        </z>
                    </fieldOfView_mm>
                </reconSpace>
                <encodingLimits>
                    <kspace_encoding_step_1>
                        <minimum>0</minimum>
                        <maximum>
                            <xsl:value-of select="siemens/YAPS/iNoOfFourierLines - 1"/>
                        </maximum>
                        <center>
<!-- Kelvin: Different logic because the shorter half of partial Fourier is acquired first -->
<!--                             <xsl:value-of select="floor(siemens/MEAS/sKSpace/lPhaseEncodingLines div 2)"/> -->
<!--                             <xsl:value-of select="ceiling((siemens/YAPS/iNoOfFourierLines - 1) div 2 * $partialFourierPhase)"/> -->
                             <xsl:value-of select="ceiling((siemens/YAPS/iNoOfFourierLines - 1) * (1 - (0.5 div $partialFourierPhase)) )"/>
                        </center>
                    </kspace_encoding_step_1>
                    <kspace_encoding_step_2>
                        <minimum>0</minimum>
                        <xsl:choose>
                            <xsl:when test="not(siemens/YAPS/iNoOfFourierPartitions) or (siemens/YAPS/i3DFTLength = 1)">
                                <maximum>0</maximum>
                                <center>0</center>
                            </xsl:when>
                            <xsl:otherwise>
                                <maximum>
                                    <xsl:value-of select="siemens/YAPS/iNoOfFourierPartitions - 1"/>
                                </maximum>
                                <xsl:choose>
                                    <xsl:when test="siemens/MEAS/sKSpace/ucTrajectory = 1">
                                        <xsl:choose>
                                            <xsl:when test="siemens/MEAS/sPat/lAccelFact3D">
                                                <xsl:choose>
                                                    <xsl:when test="not(siemens/MEAS/sPat/lAccelFact3D) > 1">
                                                        <center>
                                                            <xsl:value-of select="floor(siemens/MEAS/sKSpace/lPartitions div 2) - (siemens/YAPS/lPartitions - siemens/YAPS/iNoOfFourierPartitions)"/>
                                                        </center>
                                                    </xsl:when>
                                                    <xsl:otherwise>
                                                        <xsl:choose>
                                                            <xsl:when test="(siemens/MEAS/sKSpace/lPartitions - siemens/YAPS/iNoOfFourierPartitions) > siemens/MEAS/sPat/lAccelFact3D">
                                                                <center>
                                                                    <xsl:value-of select="floor(siemens/MEAS/sKSpace/lPartitions div 2) - (siemens/MEAS/sKSpace/lPartitions - siemens/YAPS/iNoOfFourierPartitions)"/>
                                                                </center>
                                                            </xsl:when>
                                                            <xsl:otherwise>
                                                                <center>
                                                                    <xsl:value-of select="floor(siemens/MEAS/sKSpace/lPartitions div 2)"/>
                                                                </center>
                                                            </xsl:otherwise>
                                                        </xsl:choose>
                                                    </xsl:otherwise>
                                                </xsl:choose>
                                            </xsl:when>
                                            <xsl:otherwise>
                                                <center>
                                                    <xsl:value-of select="floor(siemens/MEAS/sKSpace/lPartitions div 2) - (siemens/MEAS/sKSpace/lPartitions - siemens/YAPS/iNoOfFourierPartitions)"/>
                                                </center>
                                            </xsl:otherwise>
                                        </xsl:choose>
                                    </xsl:when>
                                    <xsl:otherwise>
                                        <center>0</center>
                                    </xsl:otherwise>
                                </xsl:choose>
                            </xsl:otherwise>
                        </xsl:choose>
                    </kspace_encoding_step_2>
                    <slice>
                        <minimum>0</minimum>
                        <maximum>
                            <xsl:value-of select="siemens/MEAS/sSliceArray/lSize - 1"/>
                        </maximum>
                        <center>0</center>
                    </slice>
                    <set>
                        <minimum>0</minimum>
                        <maximum>
                            <xsl:choose>
                                <xsl:when test="siemens/YAPS/iNSet">
                                    <xsl:value-of select="siemens/YAPS/iNSet - 1"/>
                                </xsl:when>
                                <xsl:otherwise>0</xsl:otherwise>
                            </xsl:choose>
                        </maximum>
                        <center>0</center>
                    </set>
                    <phase>
                        <minimum>0</minimum>
                        <maximum>
                            <xsl:choose>
                                <xsl:when test="siemens/MEAS/sPhysioImaging/lPhases">
                                    <xsl:value-of select="siemens/MEAS/sPhysioImaging/lPhases - 1"/>
                                </xsl:when>
                                <xsl:otherwise>0</xsl:otherwise>
                            </xsl:choose>
                        </maximum>
                        <center>0</center>
                    </phase>
                    <repetition>
                        <minimum>0</minimum>
                        <maximum>
                            <xsl:choose>
                                <xsl:when test="siemens/MEAS/lRepetitions">
                                    <xsl:value-of select="siemens/MEAS/lRepetitions"/>
                                </xsl:when>
                                <xsl:otherwise>0</xsl:otherwise>
                            </xsl:choose>
                        </maximum>
                        <center>0</center>
                    </repetition>
                    <segment>
                        <minimum>0</minimum>
                        <maximum>
                            <xsl:choose>
                                <xsl:when test="siemens/MEAS/sFastImaging/ucSegmentationMode" >
                                    <xsl:choose>
                                        <xsl:when test="siemens/MEAS/sFastImaging/ucSegmentationMode = 2">
                                            <xsl:choose>
                                                <xsl:when test="siemens/MEAS/sFastImaging/lShots">
                                                    <xsl:value-of select="siemens/MEAS/sFastImaging/lShots - 1"/>
                                                </xsl:when>
                                                <xsl:otherwise>0</xsl:otherwise>
                                            </xsl:choose>
                                        </xsl:when>
                                        <xsl:when test="siemens/MEAS/sFastImaging/ucSegmentationMode = 1">
                                            <xsl:choose>
                                                <xsl:when test="siemens/MEAS/sFastImaging/lSegments &gt; 1">
                                                    <xsl:value-of select="ceiling((siemens/YAPS/iNoOfFourierPartitions * siemens/YAPS/iNoOfFourierLines) div siemens/MEAS/sFastImaging/lSegments)"/>
                                                </xsl:when>
                                                <xsl:otherwise>0</xsl:otherwise>
                                            </xsl:choose>
                                        </xsl:when>
                                        <xsl:otherwise>0</xsl:otherwise>
                                    </xsl:choose>
                                </xsl:when>
                                <xsl:otherwise>0</xsl:otherwise>
                            </xsl:choose>
                        </maximum>
                        <center>0</center>
                    </segment>
                    <contrast>
                        <minimum>0</minimum>
                        <maximum>
                            <xsl:choose>
                                <xsl:when test="siemens/MEAS/lContrasts">
                                    <xsl:value-of select="siemens/MEAS/lContrasts - 1"/>
                                </xsl:when>
                                <xsl:otherwise>0</xsl:otherwise>
                            </xsl:choose>
                        </maximum>
                        <center>0</center>
                    </contrast>
                    <average>
                        <minimum>0</minimum>
                        <maximum>
                            <xsl:choose>
                                <xsl:when test="siemens/MEAS/lAverages">
                                    <xsl:value-of select="siemens/MEAS/lAverages - 1"/>
                                </xsl:when>
                                <xsl:otherwise>0</xsl:otherwise>
                            </xsl:choose>
                        </maximum>
                        <center>0</center>
                    </average>
                </encodingLimits>
        <parallelImaging>
          <accelerationFactor>
            <!-- Kelvin: Used a fixed acceleration rate of 3 here -->
            <kspace_encoding_step_1>3</kspace_encoding_step_1>
            <kspace_encoding_step_2>
              <xsl:choose>
                <xsl:when test="not(siemens/MEAS/sPat/lAccelFact3D)">1</xsl:when>
                <xsl:otherwise>
                  <xsl:value-of select="(siemens/MEAS/sPat/lAccelFact3D)"/>
                </xsl:otherwise>
              </xsl:choose>
            </kspace_encoding_step_2>
          </accelerationFactor>
          <calibrationMode>
                    <xsl:choose>
              <xsl:when test="siemens/MEAS/sPat/ucRefScanMode = 1">other</xsl:when>
              <xsl:when test="siemens/MEAS/sPat/ucRefScanMode = 2">embedded</xsl:when>
              <xsl:when test="siemens/MEAS/sPat/ucRefScanMode = 4">separate</xsl:when>
              <xsl:when test="siemens/MEAS/sPat/ucRefScanMode = 8">separate</xsl:when>
              <xsl:when test="siemens/MEAS/sPat/ucRefScanMode = 16">interleaved</xsl:when>
              <xsl:when test="siemens/MEAS/sPat/ucRefScanMode = 32">interleaved</xsl:when>
              <xsl:when test="siemens/MEAS/sPat/ucRefScanMode = 64">interleaved</xsl:when>
              <xsl:otherwise>other</xsl:otherwise>
                    </xsl:choose>
          </calibrationMode>
          <xsl:if test="(siemens/MEAS/sPat/ucRefScanMode = 1) or (siemens/MEAS/sPat/ucRefScanMode = 16) or (siemens/MEAS/sPat/ucRefScanMode = 32) or (siemens/MEAS/sPat/ucRefScanMode = 64)">
                    <interleavingDimension>
              <xsl:choose>
            <xsl:when test="siemens/MEAS/sPat/ucRefScanMode = 16">average</xsl:when>
            <xsl:when test="siemens/MEAS/sPat/ucRefScanMode = 32">repetition</xsl:when>
            <xsl:when test="siemens/MEAS/sPat/ucRefScanMode = 64">phase</xsl:when>
            <xsl:otherwise>other</xsl:otherwise>
              </xsl:choose>
                    </interleavingDimension>
          </xsl:if>
        </parallelImaging>
            </encoding>

            <sequenceParameters>
                <xsl:for-each select="siemens/MEAS/alTR">
                     <xsl:if test="position() = 1">
                            <TR>
                                <xsl:value-of select=". div 1000.0" />
                            </TR>
                    </xsl:if>
                    <xsl:if test="(position() &gt; 1) and (. &gt; 0)">
                        <TR>
                            <xsl:value-of select=". div 1000.0" />
                        </TR>
                    </xsl:if>
                </xsl:for-each>
                <xsl:for-each select="siemens/MEAS/alTE">
                     <xsl:if test="position() = 1">
                            <TE>
                                <xsl:value-of select=". div 1000.0" />
                            </TE>
                    </xsl:if>
                    <xsl:if test="(position() &gt; 1) and (. &gt; 0)">
                        <xsl:if test="position() &lt; ($numberOfContrasts + 1)">
                            <TE>
                                <xsl:value-of select=". div 1000.0" />
                            </TE>
                        </xsl:if>
                    </xsl:if>
                </xsl:for-each>
                <xsl:for-each select="siemens/MEAS/alTI">
                    <xsl:if test=". &gt; 0">
                        <TI>
                            <xsl:value-of select=". div 1000.0" />
                        </TI>
                    </xsl:if>
                </xsl:for-each> 
                <xsl:for-each select="siemens/DICOM/adFlipAngleDegree">
                <xsl:if test=". &gt; 0">
                        <flipAngle_deg>
                            <xsl:value-of select="." />
                        </flipAngle_deg>
                    </xsl:if>
                </xsl:for-each>
                <xsl:if test="siemens/MEAS/ucSequenceType">
                    <sequence_type>
                        <xsl:choose>
                            <xsl:when test="siemens/MEAS/ucSequenceType = 1">Flash</xsl:when>
                            <xsl:when test="siemens/MEAS/ucSequenceType = 2">SSFP</xsl:when>
                            <xsl:when test="siemens/MEAS/ucSequenceType = 4">EPI</xsl:when>
                            <xsl:when test="siemens/MEAS/ucSequenceType = 8">TurboSpinEcho</xsl:when>
                            <xsl:when test="siemens/MEAS/ucSequenceType = 16">ChemicalShiftImaging</xsl:when>
                            <xsl:when test="siemens/MEAS/ucSequenceType = 32">FID</xsl:when>
                            <xsl:otherwise>Unknown</xsl:otherwise>
                        </xsl:choose>
                    </sequence_type>
                </xsl:if>
                <xsl:if test="siemens/YAPS/lEchoSpacing">
                    <echo_spacing>
                        <xsl:value-of select="siemens/YAPS/lEchoSpacing div 1000.0" />
                    </echo_spacing>
                </xsl:if>

                <xsl:if test="siemens/YAPS/lEchoSpacing">
                    <echo_spacing>
                        <xsl:value-of select="siemens/YAPS/lEchoSpacing div 1000.0" />
                    </echo_spacing>
                </xsl:if>

                <!-- Kelvin's exported BEAT echo spacing -->
                <xsl:if test="siemens/MEAS/sWipMemBlock/alFree[63]">
                    <echo_spacing>
                      <xsl:value-of select="siemens/MEAS/sWipMemBlock/alFree[63]" />
                    </echo_spacing>
                </xsl:if>

            </sequenceParameters>

            <userParameters>
                <xsl:if test="siemens/MEAS/sAngio/sFlowArray/lSize">
                    <userParameterLong>
                        <name>VENC_0</name>
                        <value>
                            <xsl:value-of select="siemens/MEAS/sAngio/sFlowArray/asElm/s0/nVelocity" />
                        </value>
                    </userParameterLong>
                </xsl:if>

                <xsl:if test="not(siemens/MEAS/sPhysioImaging/lSignal1 = 1) and not(siemens/MEAS/sPhysioImaging/lSignal1 = 16) and (siemens/MEAS/sPhysioImaging/lMethod1 = 8)">
                    <xsl:if test="siemens/MEAS/sFastImaging/lShots >= 1">
                        <xsl:if test="siemens/MEAS/sPhysioImaging/lPhases > 1">
                            <xsl:if test="siemens/MEAS/sPhysioImaging/lRetroGatedImages > 0">
                                <userParameterLong>
                                    <name>RetroGatedImages</name>
                                    <value>
                                        <xsl:value-of select="siemens/MEAS/sPhysioImaging/lRetroGatedImages" />
                                    </value>
                                </userParameterLong>

                                <userParameterLong>
                                    <name>RetroGatedSegmentSize</name>
                                    <value>
                                        <xsl:choose>
                                            <xsl:when test="siemens/MEAS/sFastImaging/lSegments">
                                                <xsl:value-of select="siemens/MEAS/sFastImaging/lSegments"/>
                                            </xsl:when>
                                            <xsl:otherwise>0</xsl:otherwise>
                                        </xsl:choose>
                                    </value>
                                </userParameterLong>
                            </xsl:if>
                        </xsl:if>
                    </xsl:if>
                </xsl:if>

                <xsl:if test="(siemens/MEAS/ucOneSeriesForAllMeas = 2) or (siemens/MEAS/ucOneSeriesForAllMeas = 8)">
                    <userParameterLong>
                        <name>MultiSeriesForSlices</name>
                        <value>
                            <xsl:value-of select="siemens/MEAS/ucOneSeriesForAllMeas" />
                        </value>
                    </userParameterLong>
                </xsl:if>

                <xsl:if test="siemens/MEAS/sPat/lRefLinesPE">
                    <userParameterLong>
                        <name>EmbeddedRefLinesE1</name>
                        <value>
                            <xsl:value-of select="siemens/MEAS/sPat/lRefLinesPE" />
                        </value>
                    </userParameterLong>
                </xsl:if>

                <xsl:if test="siemens/MEAS/sPat/lRefLines3D">
                    <userParameterLong>
                        <name>EmbeddedRefLinesE2</name>
                        <value>
                            <xsl:value-of select="siemens/MEAS/sPat/lRefLines3D" />
                        </value>
                    </userParameterLong>
                </xsl:if>

                <xsl:if test="siemens/MEAS/lProtonDensMap">
                    <userParameterLong>
                        <name>NumOfProtonDensityImages</name>
                        <value>
                            <xsl:value-of select="siemens/MEAS/lProtonDensMap" />
                        </value>
                    </userParameterLong>
                </xsl:if>

                <xsl:choose>
                    <xsl:when test="siemens/MEAS/ucMotionCorr">
                        <userParameterLong>
                            <name>MotionCorrection</name>
                            <value>
                                <xsl:value-of select="siemens/MEAS/ucMotionCorr" />
                            </value>
                        </userParameterLong>
                    </xsl:when>
                    <xsl:otherwise>
                        <userParameterLong>
                            <name>MotionCorrection</name>
                            <value>0</value>
                        </userParameterLong>
                    </xsl:otherwise>
                </xsl:choose>

                <xsl:if test="siemens/YAPS/aflMaxwellCoefficients[1]">
                  <userParameterDouble>
                      <name>MaxwellCoefficient_0</name>
                      <value>
                          <xsl:value-of select="siemens/YAPS/aflMaxwellCoefficients[1]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>

                <xsl:if test="siemens/YAPS/aflMaxwellCoefficients[2]">
                  <userParameterDouble>
                    <name>MaxwellCoefficient_1</name>
                    <value>
                        <xsl:value-of select="siemens/YAPS/aflMaxwellCoefficients[2]" />
                    </value>
                  </userParameterDouble>
                </xsl:if>

                <xsl:if test="siemens/YAPS/aflMaxwellCoefficients[3]">
                  <userParameterDouble>
                    <name>MaxwellCoefficient_2</name>
                    <value>
                        <xsl:value-of select="siemens/YAPS/aflMaxwellCoefficients[3]" />
                    </value>
                  </userParameterDouble>
                </xsl:if>

                <xsl:if test="siemens/YAPS/aflMaxwellCoefficients[4]">
                  <userParameterDouble>
                    <name>MaxwellCoefficient_3</name>
                    <value>
                        <xsl:value-of select="siemens/YAPS/aflMaxwellCoefficients[4]" />
                    </value>
                  </userParameterDouble>
                </xsl:if>

                <!-- T2p duration (RF only) -->
                <xsl:if test="siemens/MEAS/sWipMemBlock/alFree[61]">
                  <userParameterDouble>
                      <name>T2pRfDuration</name>
                      <value>
                        <xsl:value-of select="siemens/MEAS/sWipMemBlock/alFree[61] div 1000.0" />
                      </value>
                  </userParameterDouble>
                </xsl:if>

                <!-- Time from T2p to center of k-space -->
                <xsl:if test="siemens/MEAS/sWipMemBlock/alFree[62]">
                  <userParameterDouble>
                      <name>TimeT2pToCenterKspace</name>
                      <value>
                        <xsl:value-of select="siemens/MEAS/sWipMemBlock/alFree[62] div 1000.0" />
                      </value>
                  </userParameterDouble>
                </xsl:if>

                <!--  T2 prep durations  -->
                <xsl:if test="siemens/MEAS/sPrepPulses/adT2PrepDuration[0]">
                  <userParameterDouble>
                      <name>T2PrepDuration_0</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sPrepPulses/adT2PrepDuration[0]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>
                <xsl:if test="siemens/MEAS/sPrepPulses/adT2PrepDuration[1]">
                  <userParameterDouble>
                      <name>T2PrepDuration_1</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sPrepPulses/adT2PrepDuration[1]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>
                <xsl:if test="siemens/MEAS/sPrepPulses/adT2PrepDuration[2]">
                  <userParameterDouble>
                      <name>T2PrepDuration_2</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sPrepPulses/adT2PrepDuration[2]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>
                <xsl:if test="siemens/MEAS/sPrepPulses/adT2PrepDuration[3]">
                  <userParameterDouble>
                      <name>T2PrepDuration_3</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sPrepPulses/adT2PrepDuration[3]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>
                <xsl:if test="siemens/MEAS/sPrepPulses/adT2PrepDuration[4]">
                  <userParameterDouble>
                      <name>T2PrepDuration_4</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sPrepPulses/adT2PrepDuration[4]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>
                <xsl:if test="siemens/MEAS/sPrepPulses/adT2PrepDuration[5]">
                  <userParameterDouble>
                      <name>T2PrepDuration_5</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sPrepPulses/adT2PrepDuration[5]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>
                <xsl:if test="siemens/MEAS/sPrepPulses/adT2PrepDuration[6]">
                  <userParameterDouble>
                      <name>T2PrepDuration_6</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sPrepPulses/adT2PrepDuration[6]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>
                <xsl:if test="siemens/MEAS/sPrepPulses/adT2PrepDuration[7]">
                  <userParameterDouble>
                      <name>T2PrepDuration_7</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sPrepPulses/adT2PrepDuration[7]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>
                <xsl:if test="siemens/MEAS/sPrepPulses/adT2PrepDuration[8]">
                  <userParameterDouble>
                      <name>T2PrepDuration_8</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sPrepPulses/adT2PrepDuration[8]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>
                <xsl:if test="siemens/MEAS/sPrepPulses/adT2PrepDuration[9]">
                  <userParameterDouble>
                      <name>T2PrepDuration_9</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sPrepPulses/adT2PrepDuration[9]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>
                <xsl:if test="siemens/MEAS/sPrepPulses/adT2PrepDuration[10]">
                  <userParameterDouble>
                      <name>T2PrepDuration_10</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sPrepPulses/adT2PrepDuration[10]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>
                <xsl:if test="siemens/MEAS/sPrepPulses/adT2PrepDuration[11]">
                  <userParameterDouble>
                      <name>T2PrepDuration_11</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sPrepPulses/adT2PrepDuration[11]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>
                <xsl:if test="siemens/MEAS/sPrepPulses/adT2PrepDuration[12]">
                  <userParameterDouble>
                      <name>T2PrepDuration_12</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sPrepPulses/adT2PrepDuration[12]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>
                <xsl:if test="siemens/MEAS/sPrepPulses/adT2PrepDuration[13]">
                  <userParameterDouble>
                      <name>T2PrepDuration_13</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sPrepPulses/adT2PrepDuration[13]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>
                <xsl:if test="siemens/MEAS/sPrepPulses/adT2PrepDuration[14]">
                  <userParameterDouble>
                      <name>T2PrepDuration_14</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sPrepPulses/adT2PrepDuration[14]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>
                <xsl:if test="siemens/MEAS/sPrepPulses/adT2PrepDuration[15]">
                  <userParameterDouble>
                      <name>T2PrepDuration_15</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sPrepPulses/adT2PrepDuration[15]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>
                <xsl:if test="siemens/MEAS/sPrepPulses/adT2PrepDuration[16]">
                  <userParameterDouble>
                      <name>T2PrepDuration_16</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sPrepPulses/adT2PrepDuration[16]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>

                <!--  T1rho (spin-lock) prep durations  -->
                <xsl:if test="siemens/MEAS/sWipMemBlock/adFree[1]">
                  <userParameterDouble>
                      <name>T1pPrepDuration_1</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sWipMemBlock/adFree[1]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>
                <xsl:if test="siemens/MEAS/sWipMemBlock/adFree[2]">
                  <userParameterDouble>
                      <name>T1pPrepDuration_2</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sWipMemBlock/adFree[2]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>
                <xsl:if test="siemens/MEAS/sWipMemBlock/adFree[3]">
                  <userParameterDouble>
                      <name>T1pPrepDuration_3</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sWipMemBlock/adFree[3]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>
                <xsl:if test="siemens/MEAS/sWipMemBlock/adFree[4]">
                  <userParameterDouble>
                      <name>T1pPrepDuration_4</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sWipMemBlock/adFree[4]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>
                <xsl:if test="siemens/MEAS/sWipMemBlock/adFree[5]">
                  <userParameterDouble>
                      <name>T1pPrepDuration_5</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sWipMemBlock/adFree[5]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>
                <xsl:if test="siemens/MEAS/sWipMemBlock/adFree[6]">
                  <userParameterDouble>
                      <name>T1pPrepDuration_6</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sWipMemBlock/adFree[6]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>
                <xsl:if test="siemens/MEAS/sWipMemBlock/adFree[7]">
                  <userParameterDouble>
                      <name>T1pPrepDuration_7</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sWipMemBlock/adFree[7]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>
                <xsl:if test="siemens/MEAS/sWipMemBlock/adFree[8]">
                  <userParameterDouble>
                      <name>T1pPrepDuration_8</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sWipMemBlock/adFree[8]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>
                <xsl:if test="siemens/MEAS/sWipMemBlock/adFree[9]">
                  <userParameterDouble>
                      <name>T1pPrepDuration_9</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sWipMemBlock/adFree[9]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>
                <xsl:if test="siemens/MEAS/sWipMemBlock/adFree[10]">
                  <userParameterDouble>
                      <name>T1pPrepDuration_10</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sWipMemBlock/adFree[10]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>
                <xsl:if test="siemens/MEAS/sWipMemBlock/adFree[11]">
                  <userParameterDouble>
                      <name>T1pPrepDuration_11</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sWipMemBlock/adFree[11]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>
                <xsl:if test="siemens/MEAS/sWipMemBlock/adFree[12]">
                  <userParameterDouble>
                      <name>T1pPrepDuration_12</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sWipMemBlock/adFree[12]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>
                <xsl:if test="siemens/MEAS/sWipMemBlock/adFree[13]">
                  <userParameterDouble>
                      <name>T1pPrepDuration_13</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sWipMemBlock/adFree[13]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>
                <xsl:if test="siemens/MEAS/sWipMemBlock/adFree[14]">
                  <userParameterDouble>
                      <name>T1pPrepDuration_14</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sWipMemBlock/adFree[14]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>
                <xsl:if test="siemens/MEAS/sWipMemBlock/adFree[15]">
                  <userParameterDouble>
                      <name>T1pPrepDuration_15</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sWipMemBlock/adFree[15]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>
                <xsl:if test="siemens/MEAS/sWipMemBlock/adFree[16]">
                  <userParameterDouble>
                      <name>T1pPrepDuration_16</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sWipMemBlock/adFree[16]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>

                <!-- T2/T1p prep duration (RF only) -->
                <xsl:if test="siemens/MEAS/sWipMemBlock/alFree[31]">
                  <userParameterDouble>
                      <name>T2pRfDuration_1</name>
                      <value>
                        <xsl:value-of select="siemens/MEAS/sWipMemBlock/alFree[31] div 1000.0" />
                      </value>
                  </userParameterDouble>
                </xsl:if>

                <xsl:if test="siemens/MEAS/sWipMemBlock/alFree[32]">
                  <userParameterDouble>
                      <name>T2pRfDuration_2</name>
                      <value>
                        <xsl:value-of select="siemens/MEAS/sWipMemBlock/alFree[32] div 1000.0" />
                      </value>
                  </userParameterDouble>
                </xsl:if>

                <xsl:if test="siemens/MEAS/sWipMemBlock/alFree[33]">
                  <userParameterDouble>
                      <name>T2pRfDuration_3</name>
                      <value>
                        <xsl:value-of select="siemens/MEAS/sWipMemBlock/alFree[33] div 1000.0" />
                      </value>
                  </userParameterDouble>
                </xsl:if>

                <xsl:if test="siemens/MEAS/sWipMemBlock/alFree[34]">
                  <userParameterDouble>
                      <name>T2pRfDuration_4</name>
                      <value>
                        <xsl:value-of select="siemens/MEAS/sWipMemBlock/alFree[34] div 1000.0" />
                      </value>
                  </userParameterDouble>
                </xsl:if>

                <xsl:if test="siemens/MEAS/sWipMemBlock/alFree[35]">
                  <userParameterDouble>
                      <name>T2pRfDuration_5</name>
                      <value>
                        <xsl:value-of select="siemens/MEAS/sWipMemBlock/alFree[35] div 1000.0" />
                      </value>
                  </userParameterDouble>
                </xsl:if>

                <xsl:if test="siemens/MEAS/sWipMemBlock/alFree[36]">
                  <userParameterDouble>
                      <name>T2pRfDuration_6</name>
                      <value>
                        <xsl:value-of select="siemens/MEAS/sWipMemBlock/alFree[36] div 1000.0" />
                      </value>
                  </userParameterDouble>
                </xsl:if>

                <xsl:if test="siemens/MEAS/sWipMemBlock/alFree[37]">
                  <userParameterDouble>
                      <name>T2pRfDuration_7</name>
                      <value>
                        <xsl:value-of select="siemens/MEAS/sWipMemBlock/alFree[37] div 1000.0" />
                      </value>
                  </userParameterDouble>
                </xsl:if>

                <xsl:if test="siemens/MEAS/sWipMemBlock/alFree[38]">
                  <userParameterDouble>
                      <name>T2pRfDuration_8</name>
                      <value>
                        <xsl:value-of select="siemens/MEAS/sWipMemBlock/alFree[38] div 1000.0" />
                      </value>
                  </userParameterDouble>
                </xsl:if>

                <xsl:if test="siemens/MEAS/sWipMemBlock/alFree[39]">
                  <userParameterDouble>
                      <name>T2pRfDuration_9</name>
                      <value>
                        <xsl:value-of select="siemens/MEAS/sWipMemBlock/alFree[39] div 1000.0" />
                      </value>
                  </userParameterDouble>
                </xsl:if>

                <xsl:if test="siemens/MEAS/sWipMemBlock/alFree[40]">
                  <userParameterDouble>
                      <name>T2pRfDuration_10</name>
                      <value>
                        <xsl:value-of select="siemens/MEAS/sWipMemBlock/alFree[40] div 1000.0" />
                      </value>
                  </userParameterDouble>
                </xsl:if>

                <xsl:if test="siemens/MEAS/sWipMemBlock/alFree[41]">
                  <userParameterDouble>
                      <name>T2pRfDuration_11</name>
                      <value>
                        <xsl:value-of select="siemens/MEAS/sWipMemBlock/alFree[41] div 1000.0" />
                      </value>
                  </userParameterDouble>
                </xsl:if>

                <xsl:if test="siemens/MEAS/sWipMemBlock/alFree[42]">
                  <userParameterDouble>
                      <name>T2pRfDuration_12</name>
                      <value>
                        <xsl:value-of select="siemens/MEAS/sWipMemBlock/alFree[42] div 1000.0" />
                      </value>
                  </userParameterDouble>
                </xsl:if>

                <xsl:if test="siemens/MEAS/sWipMemBlock/alFree[43]">
                  <userParameterDouble>
                      <name>T2pRfDuration_13</name>
                      <value>
                        <xsl:value-of select="siemens/MEAS/sWipMemBlock/alFree[43] div 1000.0" />
                      </value>
                  </userParameterDouble>
                </xsl:if>

                <xsl:if test="siemens/MEAS/sWipMemBlock/alFree[44]">
                  <userParameterDouble>
                      <name>T2pRfDuration_14</name>
                      <value>
                        <xsl:value-of select="siemens/MEAS/sWipMemBlock/alFree[44] div 1000.0" />
                      </value>
                  </userParameterDouble>
                </xsl:if>

                <xsl:if test="siemens/MEAS/sWipMemBlock/alFree[45]">
                  <userParameterDouble>
                      <name>T2pRfDuration_15</name>
                      <value>
                        <xsl:value-of select="siemens/MEAS/sWipMemBlock/alFree[45] div 1000.0" />
                      </value>
                  </userParameterDouble>
                </xsl:if>

                <xsl:if test="siemens/MEAS/sWipMemBlock/alFree[46]">
                  <userParameterDouble>
                      <name>T2pRfDuration_16</name>
                      <value>
                        <xsl:value-of select="siemens/MEAS/sWipMemBlock/alFree[46] div 1000.0" />
                      </value>
                  </userParameterDouble>
                </xsl:if>

                <!-- Saturation recovery times -->
                <xsl:if test="siemens/MEAS/sWipMemBlock/alFree[15]">
                  <userParameterDouble>
                      <name>SatRecTime_1</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sWipMemBlock/alFree[15]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>

                <xsl:if test="siemens/MEAS/sWipMemBlock/alFree[16]">
                  <userParameterDouble>
                      <name>SatRecTime_2</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sWipMemBlock/alFree[16]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>

                <xsl:if test="siemens/MEAS/sWipMemBlock/alFree[17]">
                  <userParameterDouble>
                      <name>SatRecTime_3</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sWipMemBlock/alFree[17]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>

                <xsl:if test="siemens/MEAS/sWipMemBlock/alFree[18]">
                  <userParameterDouble>
                      <name>SatRecTime_4</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sWipMemBlock/alFree[18]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>

                <xsl:if test="siemens/MEAS/sWipMemBlock/alFree[19]">
                  <userParameterDouble>
                      <name>SatRecTime_5</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sWipMemBlock/alFree[19]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>

                <xsl:if test="siemens/MEAS/sWipMemBlock/alFree[20]">
                  <userParameterDouble>
                      <name>SatRecTime_6</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sWipMemBlock/alFree[20]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>

                <xsl:if test="siemens/MEAS/sWipMemBlock/alFree[21]">
                  <userParameterDouble>
                      <name>SatRecTime_7</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sWipMemBlock/alFree[21]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>

                <xsl:if test="siemens/MEAS/sWipMemBlock/alFree[22]">
                  <userParameterDouble>
                      <name>SatRecTime_8</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sWipMemBlock/alFree[22]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>

                <xsl:if test="siemens/MEAS/sWipMemBlock/alFree[23]">
                  <userParameterDouble>
                      <name>SatRecTime_9</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sWipMemBlock/alFree[23]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>

                <xsl:if test="siemens/MEAS/sWipMemBlock/alFree[24]">
                  <userParameterDouble>
                      <name>SatRecTime_10</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sWipMemBlock/alFree[24]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>

                <xsl:if test="siemens/MEAS/sWipMemBlock/alFree[25]">
                  <userParameterDouble>
                      <name>SatRecTime_11</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sWipMemBlock/alFree[25]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>

                <xsl:if test="siemens/MEAS/sWipMemBlock/alFree[26]">
                  <userParameterDouble>
                      <name>SatRecTime_12</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sWipMemBlock/alFree[26]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>

                <xsl:if test="siemens/MEAS/sWipMemBlock/alFree[27]">
                  <userParameterDouble>
                      <name>SatRecTime_13</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sWipMemBlock/alFree[27]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>

                <xsl:if test="siemens/MEAS/sWipMemBlock/alFree[28]">
                  <userParameterDouble>
                      <name>SatRecTime_14</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sWipMemBlock/alFree[28]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>

                <xsl:if test="siemens/MEAS/sWipMemBlock/alFree[29]">
                  <userParameterDouble>
                      <name>SatRecTime_15</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sWipMemBlock/alFree[29]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>

                <xsl:if test="siemens/MEAS/sWipMemBlock/alFree[30]">
                  <userParameterDouble>
                      <name>SatRecTime_16</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sWipMemBlock/alFree[30]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>

                <xsl:if test="siemens/MEAS/sWipMemBlock/adFree[3]">
                  <userParameterDouble>
                      <name>HCT</name>
                      <value>
                          <xsl:value-of select="siemens/MEAS/sWipMemBlock/adFree[3]" />
                      </value>
                  </userParameterDouble>
                </xsl:if>
                
                <xsl:if test="siemens/YAPS/flContrastBolusVolume">
                    <userParameterDouble>
                        <name>ContrastBolusVolume</name>
                        <value>
                            <xsl:value-of select="siemens/YAPS/flContrastBolusVolume" />
                        </value>
                    </userParameterDouble>
                </xsl:if>

                <xsl:if test="siemens/YAPS/flContrastBolusTotalDose">
                    <userParameterDouble>
                        <name>ContrastBolusTotalDose</name>
                        <value>
                            <xsl:value-of select="siemens/YAPS/flContrastBolusTotalDose" />
                        </value>
                    </userParameterDouble>
                </xsl:if>

                <xsl:if test="siemens/YAPS/flUsedPatientWeight > 0">
                    <userParameterDouble>
                        <name>PatientWeight</name>
                        <value>
                            <xsl:value-of select="siemens/YAPS/flUsedPatientWeight" />
                        </value>
                    </userParameterDouble>
                </xsl:if>

                <xsl:if test="siemens/YAPS/flPatientHeight">
                    <userParameterDouble>
                        <name>PatientHeight</name>
                        <value>
                            <xsl:value-of select="siemens/YAPS/flPatientHeight" />
                        </value>
                    </userParameterDouble>
                </xsl:if>

                <xsl:if test="siemens/YAPS/flPatientAge">
                    <xsl:choose>
                        <xsl:when test="siemens/YAPS/flPatientAge &lt; 89">
                            <userParameterDouble>
                                <name>PatientAge</name>
                                <value>
                                    <xsl:value-of select="siemens/YAPS/flPatientAge" />
                                </value>
                            </userParameterDouble>
                        </xsl:when>
                        <xsl:otherwise>
                            <userParameterDouble>
                                <name>PatientAge</name>
                                <value>89</value>
                            </userParameterDouble>
                        </xsl:otherwise>
                    </xsl:choose>
                </xsl:if>

                <xsl:if test="siemens/YAPS/tContrastBolusAgent">
                    <userParameterString>
                        <name>ContrastBolusAgent</name>
                        <value>
                            <xsl:value-of select="siemens/YAPS/tContrastBolusAgent" />
                        </value>
                    </userParameterString>
                </xsl:if>

            </userParameters>

        </ismrmrdHeader>
    </xsl:template>

</xsl:stylesheet>
