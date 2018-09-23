#############################################################################
##  
##  PackageInfo.g for the package `PL'
##  Для механизма LoadingPackage в GAP >=4.5 минимальные установки которые нужны
##  это .PackageName, .Version и .AvailabilityTest. И если любой из этих пунктов
##  отсутствует будет случаться ошибка. Другие важные пункты это .PackageDoc и
##  .Dependencies. Остальные пункты уместны если пакет будет предоставлен для
##  других пользователей GAP, в частности если он поставляем через вэб-сайт
##  GAP'а.

SetPackageInfo( rec(

PackageName := "PL",
Subtitle := "PL: We gonna eat your brain.",
Version := "2.8.0",
Date := "29/05/2018",
PackageWWWHome :=
  Concatenation( "http://sourceforge.net/projects/plgap",
      LowercaseString( ~.PackageName ), "/" ),
ArchiveURL := Concatenation( ~.PackageWWWHome, "PL-", ~.Version ),
ArchiveFormats := ".tar.gz",

Persons := [
  rec( 
    LastName      := "Korepanov",
    FirstNames    := "Igor",
    IsAuthor      := true,
    IsMaintainer  := false,
    Email         := "___korepanov___@mail.ru",
    WWWHome       := "http://",
    PostalAddress := Concatenation( [
                       "Московский Государственный Университет Приборостроения и Информатики\n",
					   "г.Москва, ул. Стромынка 20, инд: 107996\n" ]),
    Place         := "Москва",
    Institution   := "МГУПИ"
  ),
  rec( 
    LastName      := "Korepanov",
    FirstNames    := "Alexey",
    IsAuthor      := true,
    IsMaintainer  := false,
    Email         := "",
    WWWHome       := "http://",
    PostalAddress := Concatenation( [
                       "University of Warwick\n",
                       "GB" ] ),
    Place         := "xxx",
    Institution   := "yyy"
  ),
  rec( 
    LastName      := "Sadykov",
    FirstNames    := "Nurlan",
    IsAuthor      := true,
    IsMaintainer  := true,
    Email         := "nonZeroDeterminant@gmail.com",
    WWWHome       := "http://",
    PostalAddress := Concatenation( [
                     "Московский Государственный Университет Приборостроения и информатики\n",
					 "г.Москва, ул.Стромынка 20, инд: 107996\n" ] ),
    Place         := "Москва",
    Institution   := "МГУПИ"
     ),  
],
##  Status information. Currently the following cases are recognized:
##    "accepted"      for successfully refereed packages
##    "submitted"     for packages submitted for the refereeing
##    "deposited"     for packages for which the GAP developers agreed 
##                    to distribute them with the core GAP system
##    "dev"           for development versions of packages 
##    "other"         for all other packages
##
Status := "dev",

README_URL := 
  Concatenation( ~.PackageWWWHome, "README" ),
PackageInfoURL := 
  Concatenation( ~.PackageWWWHome, "PackageInfo.g" ),
AbstractHTML := 
  "The <span class=\"pkgname\">PL</span> package, provides a collection of \
  functions for computing on PL-manifolds.",

##  Here is the information on the help books of the package, used for
##  loading into GAP's online help and maybe for an online copy of the 
##  documentation on the GAP website.
##  
##  For the online help the following is needed:
##       - the name of the book (.BookName)
##       - a long title, shown by ?books (.LongTitle, optional)
##       - the path to the manual.six file for this book (.SixFile)
##  
##  For an online version on a Web page further entries are needed,
##  if possible, provide an HTML- and a PDF-version:
##      - if there is an HTML-version the path to the start file,
##        relative to the package home directory (.HTMLStart)
##      - if there is a PDF-version the path to the .pdf-file,
##        relative to the package home directory (.PDFFile)
##      - give the paths to the files inside your package directory
##        which are needed for the online manual (as a list 
##        .ArchiveURLSubset of names of directories and/or files which 
##        should be copied from your package archive, given in .ArchiveURL 
##        above (in most cases, ["doc"] or ["doc","htm"] suffices).
##  
##  For links to other GAP or package manuals you can assume a relative 
##  position of the files as in a standard GAP installation.

##  Для ссылок на другие руководства GAP-пакетов вы можете считать относительное
##  положение этих файлов как стандартную GAP установку.

##  
# in case of several help books give a list of such records here:
# в случае нескольких справочних книг дайсте списко таких записей здесь
PackageDoc := rec(
  # use same as in GAP
  BookName  := "PL",
  # format/extension can be one of .tar.gz, .tar.bz2, -win.zip, .zoo.
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0.html",
  PDFFile   := "doc/manPL.pdf",
  # the path to the .six file used by GAP's help system
  SixFile   := "doc/manual.six",
  # a longer title of the book, this together with the book name should
  # fit on a single text line (appears with the '?books' command in GAP)
  LongTitle := "PL-manifold of a GAP Package",
),


##  Are there restrictions on the operating system for this package? Or does
##  the package need other packages to be available?
##  Существуют ограничения на операционную систему для вашего пакета? Или может
##  быть пакету нужен другой доступный пакет?
Dependencies := rec(
  # GAP version, use the version string for specifying a least version,
  # prepend a '=' for specifying an exact version.
  GAP := "4.5.3",

  # list of pairs [package name, version], package name is case
  # insensitive, exact version denoted with '=' prepended to version string.
  # without these, the package will not load
  # NeededOtherPackages := [["GAPDoc", "1.5"]],
  NeededOtherPackages := [["GAPDoc", "1.5"]],

  # list of pairs [package name, version] as above,
  # these package are will be loaded if they are available,
  # but the current package will be loaded if they are not available
  # SuggestedOtherPackages := [],
  SuggestedOtherPackages := [],

  # *Optional*: a list of pairs as above, denoting those needed packages
  # that must be completely loaded before loading of the current package
  # is started (if this is not possible due to a cyclic dependency
  # then the current package is regarded as not loadable);
  # this component should be used only if functions from the needed packages
  # in question are called (or global lists or records are accessed)
  # while the current package gets loaded
  # OtherPackagesLoadedInAdvance := [],

  # needed external conditions (programs, operating system, ...)  provide 
  # just strings as text or
  # pairs [text, URL] where URL  provides further information
  # about that point.
  # (no automatic test will be done for this, do this in your 
  # 'AvailabilityTest' function below)
  # ExternalConditions := []
  ExternalConditions := []
                      
),

##  Предоставляет тестирующую функцию для проверки готовности этого пакета. Для
##  пакетов которые содержат только GAP код, просто скажите 'ReturnTrue' здесь. 
##  Для пакетов которые могут не работать или будут иметь только часть
##  функционала, используйте 'LogPackageLodingMessange( PACKAGE_WARNING, ... )'
##  делающей запись истории сообщений которые могут быть рассмотрены позже с
##  'DisplayPackageLoadingLog'. Не вызывайте 'Print' или 'Info' в функции
##  'AvailabilityTest' пакета.

##
##  With the package loading mechanism of GAP >=4.4, the availability
##  tests of other packages, as given under .Dependencies above, will be 
##  done automatically and need not be included in this function.
##
AvailabilityTest := ReturnTrue,
#AvailabilityTest := function()
#  local path, file;
#    # test for existence of the compiled binary
#    path:= DirectoriesPackagePrograms( "example" );
#    file:= Filename( path, "hello" );
#    if file = fail then
#      LogPackageLoadingMessage( PACKAGE_WARNING,
#          [ "The program `hello' is not compiled,",
#            "`HelloWorld()' is thus unavailable.",
#            "See the installation instructions;",
#            "type: ?Installing the Example package" ] );
#    fi;
#    # if the hello binary was vital to the package we would return
#    # the following ...
#    # return file <> fail;
#    # since the hello binary is not vital we return ...
#    return true;
#  end,



##  *Optional*: path relative to package root to a file which 
##  shall be read immediately before the package is loaded.
#PreloadFile := "...",

##  *Optional*: the LoadPackage mechanism can produce a default banner from
##  the info in this file. If you are not happy with it, you can provide
##  a string here that is used as a banner. GAP decides when the banner is 
##  shown and when it is not shown (note the ~-syntax in this example).
BannerString := Concatenation( 
    "----------------------------------------------------------------\n",
    "Loading  PL ", ~.Version, "\n",
    "by ",
    JoinStringsWithSeparator( List( Filtered( ~.Persons, r -> r.IsAuthor ),
                                    r -> Concatenation(
        r.FirstNames, " ", r.LastName, " (", r.WWWHome, ")\n" ) ), "   " ),
    # "При поддержке гранта РФФИ мол_а № 14-01-31019\n",
    "----------------------------------------------------------------\n" ),

##  *Optional*, but recommended: path relative to package root to a file which 
##  contains as many tests of the package functionality as sensible.
##  The file can either consist of 'ReadTest' calls or it is itself read via
##  'ReadTest'; it is assumed that the latter case occurs if and only if
##  the file contains the string 'gap> START_TEST('.
##  For deposited packages, these tests are run regularly, as a part of the
##  standard GAP test suite.
TestFile := "tst/testall.tst",

##  *Optional*: Here you can list some keyword related to the topic 
##  of the package.
# Keywords := ["Smith normal form", "p-adic", "rational matrix inversion"]
Keywords := ["PL manifolds", "piecewise complex", "Grassmann variables",
			"ball complexes"]

));

