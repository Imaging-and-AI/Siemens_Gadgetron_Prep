<?xml version="1.0" encoding="utf-8"?>

<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:output method="xml" indent="yes"/>

  <!-- Older data with a single RF duration -->
  <xsl:variable name="singleRFDuration"
    select="contains(string(siemens/MEAS/tSequenceFileName), 'BEAT_map_PK')
         or contains(string(siemens/MEAS/tSequenceFileName), 'BEAT_map_1041C_PK2')
         or contains(string(siemens/MEAS/tSequenceFileName), 'BEAT_map_1041C_T1T2')" />

  <xsl:variable name="numarisVersion">
    <xsl:choose>
      <xsl:when test="substring(siemens/DICOM/SoftwareVersions, 0, 12) = 'syngo MR XA'">NX</xsl:when>
      <xsl:otherwise>N4</xsl:otherwise>
    </xsl:choose>
  </xsl:variable>

  <xsl:variable name="phaseOversampling">
    <xsl:choose>
      <xsl:when test="siemens/IRIS/DERIVED/phaseOversampling">
        <xsl:choose>
          <xsl:when test="string(number(siemens/IRIS/DERIVED/phaseOversampling)) = 'NaN'">0</xsl:when>
          <xsl:otherwise>
            <xsl:value-of select="siemens/IRIS/DERIVED/phaseOversampling"/>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:when>
      <!-- IRIS derived may be computed on the fly and not available, so fall back to DICOM value -->
      <xsl:when test="siemens/DICOM/PhaseOversampling">
        <xsl:choose>
          <xsl:when test="string(number(siemens/DICOM/PhaseOversampling)) = 'NaN'">0</xsl:when>
          <xsl:otherwise>
            <xsl:value-of select="siemens/DICOM/PhaseOversampling"/>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:when>
      <xsl:otherwise>0</xsl:otherwise>
    </xsl:choose>
  </xsl:variable>

  <xsl:variable name="sliceOversampling">
    <xsl:choose>
      <xsl:when test="siemens/MEAS/sKSpace/dSliceOversamplingForDialog">
        <xsl:choose>
          <xsl:when test="string(number(siemens/MEAS/sKSpace/dSliceOversamplingForDialog)) = 'NaN'">0</xsl:when>
          <xsl:otherwise>
            <xsl:value-of select="siemens/MEAS/sKSpace/dSliceOversamplingForDialog"/>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:when>
      <xsl:otherwise>0</xsl:otherwise>
    </xsl:choose>
  </xsl:variable>

  <xsl:variable name="partialFourierPhase">
    <xsl:choose>
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
    <xsl:choose>
      <xsl:when test="$numarisVersion = 'NX'">
        <xsl:value-of select="substring(siemens/IRIS/RECOMPOSE/StudyLOID, 29)"/>
      </xsl:when>
      <xsl:when test="$numarisVersion = 'N4'">
        <xsl:value-of select="substring(siemens/IRIS/RECOMPOSE/StudyLOID, 6)"/>
      </xsl:when>
    </xsl:choose>
  </xsl:variable>

  <xsl:variable name="patientID">
    <xsl:choose>
      <xsl:when test="$numarisVersion = 'NX'">
        <xsl:value-of select="substring(siemens/IRIS/RECOMPOSE/StudyLOID, 29)"/>
      </xsl:when>
      <xsl:when test="$numarisVersion = 'N4'">
        <xsl:value-of select="substring(siemens/IRIS/RECOMPOSE/PatientLOID, 6)"/>
      </xsl:when>
    </xsl:choose>
  </xsl:variable>

  <xsl:variable name="strSeperator">_</xsl:variable>

  <xsl:variable name="pixelSpacing">
    <xsl:value-of select="siemens/MEAS/sSliceArray/asSlice/s0/dReadoutFOV * 0.5 * siemens/YAPS/flReadoutOSFactor div siemens/MEAS/sKSpace/lBaseResolution"/>
  </xsl:variable>

  <xsl:template match="/">
    <ismrmrdHeader xsi:schemaLocation="http://www.ismrm.org/ISMRMRD ismrmrd.xsd"
                   xmlns="http://www.ismrm.org/ISMRMRD"
                   xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
                   xmlns:xs="http://www.w3.org/2001/XMLSchema">

      <subjectInformation>
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

        <sequenceName>
          <xsl:value-of select="siemens/MEAS/tSequenceFileName"/>
        </sequenceName>

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
              <xsl:sort data-type="number" select="." />
              <xsl:variable name="CurADC" select="."/>
              <xsl:variable name="CurADCIndex" select="position()" />
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
              <xsl:variable name="CurCoil" select="position()"/>
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
                <xsl:choose>
                  <xsl:when test="(siemens/IRIS/DERIVED/imageColumns) and (siemens/IRIS/DERIVED/imageColumns > 0)">
                    <x>
                      <xsl:value-of select="siemens/IRIS/DERIVED/imageColumns"/>
                    </x>
                  </xsl:when>
                  <xsl:otherwise>
                    <x>
                      <xsl:value-of select="siemens/MEAS/sKSpace/lBaseResolution"/>
                    </x>
                  </xsl:otherwise>
                </xsl:choose>
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
            <xsl:choose>
              <xsl:when test="(siemens/IRIS/DERIVED/imageColumns) and (siemens/IRIS/DERIVED/imageColumns > 0)">
                <x>
                  <xsl:value-of select="siemens/IRIS/DERIVED/imageColumns"/>
                </x>
              </xsl:when>
              <xsl:otherwise>
                <x>
                  <xsl:value-of select="siemens/MEAS/sKSpace/lBaseResolution"/>
                </x>
              </xsl:otherwise>
            </xsl:choose>

            <xsl:choose>
              <xsl:when test="(siemens/IRIS/DERIVED/imageLines) and (siemens/IRIS/DERIVED/imageLines > 0)">
                <y>
                  <xsl:value-of select="siemens/IRIS/DERIVED/imageLines"/>
                </y>
              </xsl:when>
              <xsl:otherwise>
                <y>
                  <xsl:value-of select="floor(siemens/MEAS/sSliceArray/asSlice/s0/dPhaseFOV div $pixelSpacing + 0.5)"/>
                </y>
              </xsl:otherwise>
            </xsl:choose>

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
              <!-- <xsl:value-of select="floor(siemens/MEAS/sKSpace/lPhaseEncodingLines div 2)"/> -->
              <!-- <xsl:value-of select="ceiling((siemens/YAPS/iNoOfFourierLines - 1) div 2 * $partialFourierPhase)"/> -->
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
                <xsl:when test="siemens/MEAS/sSpecPara/lVectorSize and siemens/MEAS/sSpecPara/lVectorSize &gt; 0">
                  <xsl:value-of select="siemens/MEAS/sKSpace/lBaseResolution - 1"/>
                </xsl:when>
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
            <center>
              <xsl:choose>
                <xsl:when test="siemens/MEAS/sSpecPara/lVectorSize and siemens/MEAS/sSpecPara/lVectorSize &gt; 0">
                  <xsl:value-of select="floor(siemens/MEAS/sKSpace/lBaseResolution div 2)"/>
                </xsl:when>
                <xsl:otherwise>0</xsl:otherwise>
              </xsl:choose>
            </center>
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

      <!-- Second encoding space: High-contrast SASHA lines acquired with acceleration rate 3 -->
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
                <xsl:choose>
                  <xsl:when test="(siemens/IRIS/DERIVED/imageColumns) and (siemens/IRIS/DERIVED/imageColumns > 0)">
                    <x>
                      <xsl:value-of select="siemens/IRIS/DERIVED/imageColumns"/>
                    </x>
                  </xsl:when>
                  <xsl:otherwise>
                    <x>
                      <xsl:value-of select="siemens/MEAS/sKSpace/lBaseResolution"/>
                    </x>
                  </xsl:otherwise>
                </xsl:choose>
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
            <xsl:choose>
              <xsl:when test="(siemens/IRIS/DERIVED/imageColumns) and (siemens/IRIS/DERIVED/imageColumns > 0)">
                <x>
                  <xsl:value-of select="siemens/IRIS/DERIVED/imageColumns"/>
                </x>
              </xsl:when>
              <xsl:otherwise>
                <x>
                  <xsl:value-of select="siemens/MEAS/sKSpace/lBaseResolution"/>
                </x>
              </xsl:otherwise>
            </xsl:choose>

            <xsl:choose>
              <xsl:when test="(siemens/IRIS/DERIVED/imageLines) and (siemens/IRIS/DERIVED/imageLines > 0)">
                <y>
                  <xsl:value-of select="siemens/IRIS/DERIVED/imageLines"/>
                </y>
              </xsl:when>
              <xsl:otherwise>
                <y>
                  <xsl:value-of select="floor(siemens/MEAS/sSliceArray/asSlice/s0/dPhaseFOV div $pixelSpacing + 0.5)"/>
                </y>
              </xsl:otherwise>
            </xsl:choose>

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
              <!-- <xsl:value-of select="floor(siemens/MEAS/sKSpace/lPhaseEncodingLines div 2)"/> -->
              <!-- <xsl:value-of select="ceiling((siemens/YAPS/iNoOfFourierLines - 1) div 2 * $partialFourierPhase)"/> -->
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
                <xsl:when test="siemens/MEAS/sSpecPara/lVectorSize and siemens/MEAS/sSpecPara/lVectorSize &gt; 0">
                  <xsl:value-of select="siemens/MEAS/sKSpace/lBaseResolution - 1"/>
                </xsl:when>
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
            <center>
              <xsl:choose>
                <xsl:when test="siemens/MEAS/sSpecPara/lVectorSize and siemens/MEAS/sSpecPara/lVectorSize &gt; 0">
                  <xsl:value-of select="floor(siemens/MEAS/sKSpace/lBaseResolution div 2)"/>
                </xsl:when>
                <xsl:otherwise>0</xsl:otherwise>
              </xsl:choose>
            </center>
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
        <xsl:if test="siemens/MEAS/ulVersion">
          <userParameterLong>
            <name>ulVersion</name>
            <value>
              <xsl:value-of select="siemens/MEAS/ulVersion" />
            </value>
          </userParameterLong>
        </xsl:if>

        <xsl:if test="siemens/MEAS/sAngio/sFlowArray/lSize">
          <xsl:if test="not(siemens/MEAS/sAngio/sFlowArray/asElm/s0/nVelocity = 0)">
            <userParameterLong>
              <name>VENC_0</name>
              <value>
                <xsl:value-of select="siemens/MEAS/sAngio/sFlowArray/asElm/s0/nVelocity" />
              </value>
            </userParameterLong>
          </xsl:if>
        </xsl:if>

        <!-- Any signal source other than respiratory, with retrogating -->
        <xsl:if test="not(siemens/MEAS/sPhysioImaging/lSignal1 = 1) and not(siemens/MEAS/sPhysioImaging/lSignal1 = 16) and (siemens/MEAS/sPhysioImaging/lMethod1 = 8)">
          <xsl:if test="siemens/MEAS/sPhysioImaging/lRetroGatedImages > 0">
            <userParameterLong>
              <name>RetroGatedImages</name>
              <value>
                <xsl:value-of select="siemens/MEAS/sPhysioImaging/lRetroGatedImages"/>
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

        <xsl:if test="(siemens/MEAS/ucOneSeriesForAllMeas = 2) or (siemens/MEAS/ucOneSeriesForAllMeas = 8)">
          <userParameterLong>
            <name>MultiSeriesForSlices</name>
            <value>
              <xsl:value-of select="siemens/MEAS/ucOneSeriesForAllMeas" />
            </value>
          </userParameterLong>
        </xsl:if>

        <xsl:if test="(siemens/MEAS/sPat/lRefLinesPE) and not(siemens/MEAS/sPat/lRefLinesPE = 0)">
          <userParameterLong>
            <name>EmbeddedRefLinesE1</name>
            <value>
              <xsl:value-of select="siemens/MEAS/sPat/lRefLinesPE" />
            </value>
          </userParameterLong>
        </xsl:if>

        <xsl:if test="siemens/MEAS/sPat/lRefLines3D and not(siemens/MEAS/sPat/lRefLines3D = 0)">
          <userParameterLong>
            <name>EmbeddedRefLinesE2</name>
            <value>
              <xsl:value-of select="siemens/MEAS/sPat/lRefLines3D" />
            </value>
          </userParameterLong>
        </xsl:if>

        <xsl:if test="siemens/MEAS/lProtonDensMap and not(siemens/MEAS/lProtonDensMap = 0)">
          <userParameterLong>
            <name>NumOfProtonDensityImages</name>
            <value>
              <xsl:value-of select="siemens/MEAS/lProtonDensMap" />
            </value>
          </userParameterLong>
        </xsl:if>

        <xsl:choose>
          <xsl:when test="siemens/MEAS/ucMotionCorr and not(siemens/MEAS/ucMotionCorr = 0)">
            <userParameterLong>
              <name>MotionCorrection</name>
              <value><xsl:value-of select="siemens/MEAS/ucMotionCorr"/></value>
            </userParameterLong>
          </xsl:when>

          <xsl:when test="siemens/MEAS/ulMotionCorr and not(siemens/MEAS/ulMotionCorr = 0)">
            <userParameterLong>
              <name>MotionCorrection</name>
              <value><xsl:value-of select="siemens/MEAS/ulMotionCorr"/></value>
            </userParameterLong>
          </xsl:when>
        </xsl:choose>

        <xsl:if test="siemens/IRIS/DERIVED/relSliceNumber">
          <xsl:for-each select="siemens/IRIS/DERIVED/relSliceNumber">
            <xsl:if test=". &gt; -1">
              <userParameterLong>
                <name>
                  <xsl:value-of select="concat('RelativeSliceNumber_', position())"/>
                </name>
                <value>
                  <xsl:value-of select="." />
                </value>
              </userParameterLong>
            </xsl:if>
          </xsl:for-each>
        </xsl:if>

        <xsl:if test="siemens/MEAS/sSpecPara/lVectorSize and siemens/MEAS/sSpecPara/lVectorSize &gt; 0">
          <userParameterLong>
            <name>SpecVectorSize</name>
            <value>
              <xsl:value-of select="siemens/MEAS/sSpecPara/lVectorSize" />
            </value>
          </userParameterLong>
        </xsl:if>

        <xsl:if test="siemens/MEAS/sSpecPara/lVectorSize and siemens/MEAS/sSpecPara/lVectorSize &gt; 0">
          <xsl:if test="siemens/MEAS/sSpecPara/ucRemoveOversampling">
            <userParameterLong>
              <name>SpecRemoveOversampling</name>
              <value>
                <xsl:choose>
                  <xsl:when test="
                    siemens/MEAS/sSpecPara/ucRemoveOversampling = 'true' or 
                    siemens/MEAS/sSpecPara/ucRemoveOversampling = '1' or
                    siemens/MEAS/sSpecPara/ucRemoveOversampling = '0x1' ">1</xsl:when>
                  <xsl:otherwise>0</xsl:otherwise>
                </xsl:choose>
              </value>
            </userParameterLong>
          </xsl:if>
        </xsl:if>

        <xsl:if test="siemens/MEAS/sDiffusion/alAverages">
          <xsl:for-each select="siemens/MEAS/sDiffusion/alAverages">
            <xsl:if test=". &gt; 0">
              <userParameterLong>
                <name>
                  <xsl:value-of select="concat('DiffusionAverages_', position())"/>
                </name>
                <value>
                  <xsl:value-of select="." />
                </value>
              </userParameterLong>
            </xsl:if>
          </xsl:for-each>
        </xsl:if>

        <!-- Translate all non-zero sWipMemBlock/alFree parameters -->
        <xsl:if test="siemens/MEAS/sWipMemBlock/alFree">
          <xsl:for-each select="siemens/MEAS/sWipMemBlock/alFree">
            <xsl:if test="not(. = 0)">
              <userParameterLong>
                <name>
                  <xsl:value-of select="concat('sWipMemBlock_alFree_', position()-1)"/>
                </name>
                <value>
                  <xsl:value-of select="." />
                </value>
              </userParameterLong>
            </xsl:if>
          </xsl:for-each>
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

        <xsl:if test="siemens/YAPS/flReadoutOSFactor">
          <userParameterDouble>
          <name>ReadoutOSFactor</name>
          <value>
            <xsl:value-of select="siemens/YAPS/flReadoutOSFactor" />
          </value>
          </userParameterDouble>
        </xsl:if>

        <xsl:if test="siemens/MEAS/sRXSPEC/alDwellTime">
          <xsl:for-each select="siemens/MEAS/sRXSPEC/alDwellTime">
            <xsl:if test="position() = 1">
              <userParameterDouble>
                <name>
                  <xsl:value-of select="concat('DwellTime_', position()-1)"/>
                </name>
                <value>
                  <xsl:value-of select=". div 1000.0" />
                </value>
              </userParameterDouble>
            </xsl:if>
            <xsl:if test="(position() &gt; 1) and (position() &lt; ($numberOfContrasts + 1)) and (. &gt; 0)">
              <userParameterDouble>
                <name>
                  <xsl:value-of select="concat('DwellTime_', position()-1)"/>
                </name>
                <value>
                  <xsl:value-of select=". div 1000.0" />
                </value>
              </userParameterDouble>
            </xsl:if>
          </xsl:for-each>
        </xsl:if>

        <xsl:if test="siemens/MEAS/sSpecPara/lVectorSize and siemens/MEAS/sSpecPara/lVectorSize &gt; 0">
          <xsl:variable name="DwellTime" select="siemens/MEAS/sRXSPEC/alDwellTime[1]"/>
          <xsl:variable name="ReadoutOSFactor" select="siemens/YAPS/flReadoutOSFactor"/>
          <xsl:variable name="SpectralWidth" select="(1000000000 div number($DwellTime) div number($ReadoutOSFactor))"/>
          <userParameterDouble>
            <name>SpectralWidth</name>
            <value>
              <xsl:value-of select="$SpectralWidth" />
            </value>
          </userParameterDouble>

          <xsl:if test="siemens/MEAS/sSpecPara/sVoI/dThickness">
            <userParameterDouble>
            <name>SpecVoiThickness</name>
            <value>
                <xsl:value-of select="siemens/MEAS/sSpecPara/sVoI/dThickness" />
            </value>
            </userParameterDouble>
          </xsl:if>

          <xsl:if test="siemens/MEAS/sSpecPara/sVoI/dPhaseFOV">
            <userParameterDouble>
            <name>SpecVoiPhaseFOV</name>
            <value>
                <xsl:value-of select="siemens/MEAS/sSpecPara/sVoI/dPhaseFOV" />
            </value>
            </userParameterDouble>
          </xsl:if>

          <xsl:if test="siemens/MEAS/sSpecPara/sVoI/dReadoutFOV">
            <userParameterDouble>
            <name>SpecVoiReadoutFOV</name>
            <value>
                <xsl:value-of select="siemens/MEAS/sSpecPara/sVoI/dReadoutFOV" />
            </value>
            </userParameterDouble>
          </xsl:if>
        </xsl:if>

        <xsl:if test="siemens/YAPS/aflMaxwellCoefficients">
          <xsl:for-each select="siemens/YAPS/aflMaxwellCoefficients">
            <xsl:if test="not(. = 0)">
              <userParameterDouble>
                <name>
                  <xsl:value-of select="concat('MaxwellCoefficient_', position()-1)"/>
                </name>
                <value>
                  <xsl:value-of select="." />
                </value>
              </userParameterDouble>
            </xsl:if>
          </xsl:for-each>
        </xsl:if>

        <!-- Translate all non-zero sWipMemBlock/adFree parameters -->
        <xsl:if test="siemens/MEAS/sWipMemBlock/adFree">
          <xsl:for-each select="siemens/MEAS/sWipMemBlock/adFree">
            <xsl:if test="not(. = 0)">
              <userParameterDouble>
                <name>
                  <xsl:value-of select="concat('sWipMemBlock_adFree_', position()-1)"/>
                </name>
                <value>
                  <xsl:value-of select="." />
                </value>
              </userParameterDouble>
            </xsl:if>
          </xsl:for-each>
        </xsl:if>

        <!-- T2p duration (RF only) -->
        <xsl:if test="siemens/MEAS/sWipMemBlock/alFree[61] and not(siemens/MEAS/sWipMemBlock/alFree[61] = 0)">
          <userParameterDouble>
              <name>T2pRfDuration</name>
              <value>
                <xsl:value-of select="siemens/MEAS/sWipMemBlock/alFree[61] div 1000.0" />
              </value>
          </userParameterDouble>
        </xsl:if>

        <!-- Time from T2p to center of k-space -->
        <xsl:if test="siemens/MEAS/sWipMemBlock/alFree[62] and not(siemens/MEAS/sWipMemBlock/alFree[62] = 0)">
          <userParameterDouble>
              <name>TimeT2pToCenterKspace</name>
              <value>
                <xsl:value-of select="siemens/MEAS/sWipMemBlock/alFree[62] div 1000.0" />
              </value>
          </userParameterDouble>
        </xsl:if>

        <!--  T2 prep durations  -->
        <xsl:if test="siemens/MEAS/sPrepPulses/adT2PrepDuration">
          <xsl:for-each select="siemens/MEAS/sPrepPulses/adT2PrepDuration">
            <userParameterDouble>
              <name>
                <xsl:value-of select="concat('T2PrepDuration_', position())"/>
              </name>
              <value>
                <xsl:value-of select="." />
              </value>
            </userParameterDouble>
          </xsl:for-each>
        </xsl:if>

        <!--  T1rho (spin-lock) prep durations  -->
        <xsl:if test="not($singleRFDuration)">
          <xsl:if test="siemens/MEAS/sWipMemBlock/adFree">
            <xsl:for-each select="siemens/MEAS/sWipMemBlock/adFree[position() &lt;= 16]">
              <userParameterDouble>
                <name>
                  <xsl:value-of select="concat('T1pPrepDuration_', position())"/>
                </name>
                <value>
                  <xsl:value-of select="." />
                </value>
              </userParameterDouble>
            </xsl:for-each>
          </xsl:if>

          <!-- T2/T1p prep duration (RF only) -->
          <xsl:if test="siemens/MEAS/sWipMemBlock/alFree">
            <xsl:for-each select="siemens/MEAS/sWipMemBlock/alFree[position() &gt;= 31 and position() &lt;= 46]">
              <userParameterDouble>
                <name>
                  <xsl:value-of select="concat('T2pRfDuration_', position())"/>
                </name>
                <value>
                  <xsl:value-of select=". div 1000.0" />
                </value>
              </userParameterDouble>
            </xsl:for-each>
          </xsl:if>
        </xsl:if>

        <!-- Saturation recovery times -->
        <xsl:if test="siemens/MEAS/sWipMemBlock/alFree">
          <xsl:for-each select="siemens/MEAS/sWipMemBlock/alFree[position() &gt;= 15 and position() &lt;= 30]">
            <userParameterDouble>
              <name>
                <xsl:value-of select="concat('SatRecTime_', position())"/>
              </name>
              <value>
                <xsl:value-of select="." />
              </value>
            </userParameterDouble>
          </xsl:for-each>
        </xsl:if>

        <xsl:if test="siemens/MEAS/sWipMemBlock/adFree[3] and not(siemens/MEAS/sWipMemBlock/adFree[3] = 0)">
          <userParameterDouble>
            <name>HCT</name>
            <value>
              <xsl:value-of select="siemens/MEAS/sWipMemBlock/adFree[3]" />
            </value>
          </userParameterDouble>
        </xsl:if>

        <xsl:if test="siemens/YAPS/flContrastBolusVolume and not(siemens/YAPS/flContrastBolusVolume = 0)">
          <userParameterDouble>
            <name>ContrastBolusVolume</name>
            <value>
              <xsl:value-of select="siemens/YAPS/flContrastBolusVolume" />
            </value>
          </userParameterDouble>
        </xsl:if>

        <xsl:if test="siemens/YAPS/flContrastBolusTotalDose and not(siemens/YAPS/flContrastBolusTotalDose = 0)">
          <userParameterDouble>
            <name>ContrastBolusTotalDose</name>
            <value>
              <xsl:value-of select="siemens/YAPS/flContrastBolusTotalDose" />
            </value>
          </userParameterDouble>
        </xsl:if>

        <xsl:if test="siemens/DICOM/SoftwareVersions">
          <userParameterString>
            <name>SoftwareVersions</name>
            <value>
              <xsl:value-of select="siemens/DICOM/SoftwareVersions" />
            </value>
          </userParameterString>
        </xsl:if>

        <xsl:if test="siemens/IRIS/RECOMPOSE/StudyLOID">
          <userParameterString>
            <name>StudyLOID</name>
            <value>
              <xsl:value-of select="siemens/IRIS/RECOMPOSE/StudyLOID" />
            </value>
          </userParameterString>
        </xsl:if>

        <!-- Don't include for WATER_SUPPRESSION_OFF case -->
        <xsl:if test="siemens/MEAS/sPrepPulses/ucWaterSat and not(siemens/MEAS/sPrepPulses/ucWaterSat = 4)">
          <userParameterString>
            <name>WaterSaturation</name>
            <value>
              <xsl:choose>
                <xsl:when test="siemens/MEAS/sPrepPulses/ucWaterSat = 1">WATER_SATURATION</xsl:when>
                <xsl:when test="siemens/MEAS/sPrepPulses/ucWaterSat = 2">FAT_EXCITATION</xsl:when>
                <xsl:when test="siemens/MEAS/sPrepPulses/ucWaterSat = 4">WATER_SUPPRESSION_OFF</xsl:when>
                <xsl:when test="siemens/MEAS/sPrepPulses/ucWaterSat = 8">WATER_SATURATION_QUICK</xsl:when>
                <xsl:when test="siemens/MEAS/sPrepPulses/ucWaterSat = 16">WATER_SUPPRESSION_PARTIAL</xsl:when>
                <xsl:when test="siemens/MEAS/sPrepPulses/ucWaterSat = 32">WATER_SUPPRESSION_WEAK</xsl:when>
                <xsl:when test="siemens/MEAS/sPrepPulses/ucWaterSat = 64">WATER_SUPPRESSION_RF_OFF</xsl:when>
                <xsl:otherwise>UNDEFINED</xsl:otherwise>
              </xsl:choose>
            </value>
          </userParameterString>
        </xsl:if>

        <!-- Don't include for FAT_SUPPRESSION_OFF case -->
        <xsl:if test="siemens/MEAS/sPrepPulses/lFatWaterContrast and not(siemens/MEAS/sPrepPulses/lFatWaterContrast = 1)">
          <userParameterString>
            <name>FatWaterContrast</name>
            <value>
              <xsl:choose>
                <xsl:when test="siemens/MEAS/sPrepPulses/lFatWaterContrast =     1">FAT_SUPPRESSION_OFF</xsl:when>
                <xsl:when test="siemens/MEAS/sPrepPulses/lFatWaterContrast =     4">FAT_SATURATION</xsl:when>
                <xsl:when test="siemens/MEAS/sPrepPulses/lFatWaterContrast =     8">FAST_FAT_SATURATION</xsl:when>
                <xsl:when test="siemens/MEAS/sPrepPulses/lFatWaterContrast =    16">WATER_EXCITATION</xsl:when>
                <xsl:when test="siemens/MEAS/sPrepPulses/lFatWaterContrast =    32">FAST_WATER_EXCITATION</xsl:when>
                <xsl:when test="siemens/MEAS/sPrepPulses/lFatWaterContrast =    64">WATER_SATURATION</xsl:when>
                <xsl:when test="siemens/MEAS/sPrepPulses/lFatWaterContrast =   128">FAT_EXCITATION</xsl:when>
                <xsl:when test="siemens/MEAS/sPrepPulses/lFatWaterContrast =   256">FAST_WATER_SATURATION</xsl:when>
                <xsl:when test="siemens/MEAS/sPrepPulses/lFatWaterContrast =   512">PARTIAL_WATER_SATURATION</xsl:when>
                <xsl:when test="siemens/MEAS/sPrepPulses/lFatWaterContrast =  1024">WATER_SUPPRESSION_WEAK</xsl:when>
                <xsl:when test="siemens/MEAS/sPrepPulses/lFatWaterContrast =  2048">WATER_SUPPRESSION_RF_OFF</xsl:when>
                <xsl:when test="siemens/MEAS/sPrepPulses/lFatWaterContrast =  4096">SPAIR_FAT_SUPPRESSION</xsl:when>
                <xsl:when test="siemens/MEAS/sPrepPulses/lFatWaterContrast =  8192">DIXON_FAT_WATER_SEPARATION</xsl:when>
                <xsl:when test="siemens/MEAS/sPrepPulses/lFatWaterContrast = 16384">FAST_DIXON_FAT_WATER_SEPARATION</xsl:when>
                <xsl:otherwise>UNDEFINED</xsl:otherwise>
              </xsl:choose>
            </value>
          </userParameterString>
        </xsl:if>

        <!-- Don't include for FAT_SUPPRESSION_OFF case -->
        <xsl:if test="siemens/MEAS/sPrepPulses/ucFatSat and not(siemens/MEAS/sPrepPulses/ucFatSat = 4)">
          <userParameterString>
            <name>FatSat</name>
            <value>
              <xsl:choose>
                <xsl:when test="siemens/MEAS/sPrepPulses/ucFatSat =  1">FAT_SATURATION</xsl:when>
                <xsl:when test="siemens/MEAS/sPrepPulses/ucFatSat =  2">WATER_EXCITATION</xsl:when>
                <xsl:when test="siemens/MEAS/sPrepPulses/ucFatSat =  8">FAT_SATURATION_QUICK</xsl:when>
                <xsl:when test="siemens/MEAS/sPrepPulses/ucFatSat = 16">WATER_EXCITATION_FAST</xsl:when>
                <xsl:when test="siemens/MEAS/sPrepPulses/ucFatSat = 32">FAT_SUPPRESSION_OPTIMAL</xsl:when>
                <xsl:otherwise>UNDEFINED</xsl:otherwise>
              </xsl:choose>
            </value>
          </userParameterString>
        </xsl:if>

        <!-- Translate sWipMemBlock/tFree parameter -->
        <xsl:if test="siemens/MEAS/sWipMemBlock/tFree and (string(siemens/MEAS/sWipMemBlock/tFree))">
          <userParameterString>
            <name>sWipMemBlock_tFree</name>
            <value>
              <xsl:value-of select="siemens/MEAS/sWipMemBlock/tFree" />
            </value>
          </userParameterString>
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