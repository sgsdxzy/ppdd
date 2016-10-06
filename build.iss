; Script generated by the Inno Setup Script Wizard.
; SEE THE DOCUMENTATION FOR DETAILS ON CREATING INNO SETUP SCRIPT FILES!

#define MyAppName "PPDD"
#define MyAppVersion "0.99"
#define MyAppPublisher "XU Zhiyi"
#define MyAppExeName "ppdd_gui.exe"

[Setup]
; NOTE: The value of AppId uniquely identifies this application.
; Do not use the same AppId value in installers for other applications.
; (To generate a new GUID, click Tools | Generate GUID inside the IDE.)
AppId={{1E372463-6779-40CA-98B7-4F69910CB934}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
;AppVerName={#MyAppName} {#MyAppVersion}
AppPublisher={#MyAppPublisher}
DefaultDirName={pf}\{#MyAppName}
DefaultGroupName={#MyAppName}
AllowNoIcons=yes
LicenseFile=C:\Users\sgsdxzy\Sync\Programs\ppdd\LICENSE
OutputBaseFilename=ppdd
SetupIconFile=C:\Users\sgsdxzy\Sync\Programs\ppdd\ppdd.ico
Compression=lzma
SolidCompression=yes

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"; Flags: unchecked

[Files]
Source: "C:\Users\sgsdxzy\Programs\ppdd\dist\ppdd_gui.exe"; DestDir: "{app}"; Flags: ignoreversion
Source: "C:\Users\sgsdxzy\Programs\ppdd\dist\application-exit.png"; DestDir: "{app}"; Flags: ignoreversion
Source: "C:\Users\sgsdxzy\Programs\ppdd\dist\document-open.png"; DestDir: "{app}"; Flags: ignoreversion
Source: "C:\Users\sgsdxzy\Programs\ppdd\dist\document-save.png"; DestDir: "{app}"; Flags: ignoreversion
Source: "C:\Users\sgsdxzy\Programs\ppdd\dist\document-save-as.png"; DestDir: "{app}"; Flags: ignoreversion
Source: "C:\Users\sgsdxzy\Programs\ppdd\dist\help-about.png"; DestDir: "{app}"; Flags: ignoreversion
Source: "C:\Users\sgsdxzy\Programs\ppdd\dist\libiomp5md.dll"; DestDir: "{app}"; Flags: ignoreversion
Source: "C:\Users\sgsdxzy\Programs\ppdd\dist\media-playback-start.png"; DestDir: "{app}"; Flags: ignoreversion
Source: "C:\Users\sgsdxzy\Programs\ppdd\dist\media-seek-forward.png"; DestDir: "{app}"; Flags: ignoreversion
Source: "C:\Users\sgsdxzy\Programs\ppdd\dist\mkl_avx.dll"; DestDir: "{app}"; Flags: ignoreversion
Source: "C:\Users\sgsdxzy\Programs\ppdd\dist\mkl_avx2.dll"; DestDir: "{app}"; Flags: ignoreversion
Source: "C:\Users\sgsdxzy\Programs\ppdd\dist\mkl_core.dll"; DestDir: "{app}"; Flags: ignoreversion
Source: "C:\Users\sgsdxzy\Programs\ppdd\dist\mkl_intel_thread.dll"; DestDir: "{app}"; Flags: ignoreversion
Source: "C:\Users\sgsdxzy\Programs\ppdd\dist\mkl_p4.dll"; DestDir: "{app}"; Flags: ignoreversion
Source: "C:\Users\sgsdxzy\Programs\ppdd\dist\mkl_p4m.dll"; DestDir: "{app}"; Flags: ignoreversion
Source: "C:\Users\sgsdxzy\Programs\ppdd\dist\mkl_p4m3.dll"; DestDir: "{app}"; Flags: ignoreversion
Source: "C:\Users\sgsdxzy\Programs\ppdd\dist\mkl_sequential.dll"; DestDir: "{app}"; Flags: ignoreversion
Source: "C:\Users\sgsdxzy\Programs\ppdd\dist\mkl_tbb_thread.dll"; DestDir: "{app}"; Flags: ignoreversion
Source: "C:\Users\sgsdxzy\Programs\ppdd\dist\mkl_vml_avx.dll"; DestDir: "{app}"; Flags: ignoreversion
Source: "C:\Users\sgsdxzy\Programs\ppdd\dist\mkl_vml_avx2.dll"; DestDir: "{app}"; Flags: ignoreversion
Source: "C:\Users\sgsdxzy\Programs\ppdd\dist\mkl_vml_cmpt.dll"; DestDir: "{app}"; Flags: ignoreversion
Source: "C:\Users\sgsdxzy\Programs\ppdd\dist\mkl_vml_ia.dll"; DestDir: "{app}"; Flags: ignoreversion
Source: "C:\Users\sgsdxzy\Programs\ppdd\dist\mkl_vml_p4.dll"; DestDir: "{app}"; Flags: ignoreversion
Source: "C:\Users\sgsdxzy\Programs\ppdd\dist\mkl_vml_p4m.dll"; DestDir: "{app}"; Flags: ignoreversion
Source: "C:\Users\sgsdxzy\Programs\ppdd\dist\mkl_vml_p4m2.dll"; DestDir: "{app}"; Flags: ignoreversion
Source: "C:\Users\sgsdxzy\Programs\ppdd\dist\mkl_vml_p4m3.dll"; DestDir: "{app}"; Flags: ignoreversion
Source: "C:\Users\sgsdxzy\Programs\ppdd\dist\process-stop.png"; DestDir: "{app}"; Flags: ignoreversion
Source: "C:\Users\sgsdxzy\Programs\ppdd\dist\python34.dll"; DestDir: "{app}"; Flags: ignoreversion
Source: "C:\Users\sgsdxzy\Programs\ppdd\dist\mpl-data\*"; DestDir: "{app}\mpl-data"; Flags: ignoreversion recursesubdirs createallsubdirs
Source: "C:\Users\sgsdxzy\Programs\ppdd\dist\platforms\*"; DestDir: "{app}\platforms"; Flags: ignoreversion recursesubdirs createallsubdirs
Source: "C:\Users\sgsdxzy\Programs\ppdd\dist\pylibs\*"; DestDir: "{app}\pylibs"; Flags: ignoreversion recursesubdirs createallsubdirs
; NOTE: Don't use "Flags: ignoreversion" on any shared system files

[Icons]
Name: "{group}\{#MyAppName}"; Filename: "{app}\{#MyAppExeName}"
Name: "{group}\{cm:UninstallProgram,{#MyAppName}}"; Filename: "{uninstallexe}"
Name: "{commondesktop}\{#MyAppName}"; Filename: "{app}\{#MyAppExeName}"; Tasks: desktopicon

[Run]
Filename: "{app}\{#MyAppExeName}"; Description: "{cm:LaunchProgram,{#StringChange(MyAppName, '&', '&&')}}"; Flags: nowait postinstall skipifsilent

