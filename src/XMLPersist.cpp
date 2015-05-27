#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include "XMLPersist.h"
#include "TimeCounter.h"
#include "MesmerConfig.h"
#include "plugin.h"

using namespace std;

namespace mesmer
{
  //static variables
  TiXmlDocument* XMLPersist::pDefaults=NULL;
  std::string  XMLPersist::m_defaultsfilename;
  std::string  XMLPersist::xtraerr("\nIt may not exist, it may be being used by another program,"
                 "\nor it may not consist of well-formed XML.");


  PersistPtr XMLPersist::XmlLoad(const std::string& inputfilename,
    const std::string& defaultsfilename,
    const std::string& title)
  {
    if(!defaultsfilename.empty())
      m_defaultsfilename = defaultsfilename;
    TiXmlDocument* pdoc = new TiXmlDocument();//Deleted in destructor

    bool ret;
    if(inputfilename.empty())
      ret = pdoc->LoadFile(stdin);
    else
      ret = pdoc->LoadFile(inputfilename);
    if( !ret )
    {
      cerr << "Could not load file " << inputfilename
           << xtraerr  << endl;

      delete pdoc;
      return PersistPtr(NULL);
    }

    TiXmlElement* root = pdoc->RootElement();
    if(title.size() && root->ValueStr()!=title)
    {
      cerr << inputfilename << " does not have a required root element named " << title << endl;
      delete pdoc;
      return PersistPtr(NULL);
    }
    return PersistPtr(new XMLPersist(root, pdoc));//this XMLPersist object has a document pointer
  }

  bool XMLPersist::XmlInclude(const std::string& filename)
  {
    //See documentation in declaration
    TiXmlDocument extradoc;
    if(!extradoc.LoadFile(filename))
    {
      cerr << "Could not load the secondary input file " << filename << xtraerr << endl;
      return false;
    }
    TiXmlNode* pnExtra = extradoc.RootElement();
    string extraRootStr(pnExtra->ValueStr());
    //Find child element of <mesmer> which matches the root element of the extra file
    TiXmlNode* pnMainSection = pnNode->FirstChild(extraRootStr);
    //if successful is case A
    if(!pnMainSection)
    {
      //case B
      //Look for the main section (e.g. moleculeList) which contains <molecule>s
      TiXmlNode* pnExtraChild = pnExtra->FirstChild();//e.g. <molecule>
      string extraChildStr(pnExtraChild->ValueStr());//"molecule"
      for(pnMainSection = pnNode->FirstChild(); pnMainSection; pnMainSection = pnMainSection->NextSibling() )
        // <moleculeList>, <reactionList>,...
      {
        if(pnMainSection->FirstChild(extraChildStr))//<moleculeList> has child <molecules>
          break;
      }
      if(!pnMainSection)
      {
        cerr << "Could not include the contents of " << filename << endl;
        return false;
      }
    }

    //insert each child element from extrafile
    TiXmlNode* pnMainBefore = pnMainSection->FirstChild();
    
    for(TiXmlNode* child = pnExtra->LastChild(); child; child = child->PreviousSibling() )
    {
      if(!pnMainBefore)
        pnMainBefore = pnMainSection->InsertEndChild(*child);
      else
        pnMainBefore = pnMainSection->InsertBeforeChild(pnMainBefore, *child);
      PersistPtr(new XMLPersist(pnMainBefore->ToElement()))->XmlWriteMetadata("source", filename);
    }

    return true;
  }

  XMLPersist::~XMLPersist()
  {
    delete pDocument; //doesn't matter that pDocument is usually NULL
  }

  PersistPtr XMLPersist::XmlMoveTo(const std::string& name) const
  {
    TiXmlElement* pnEl = pnNode->FirstChildElement(name);
    if(!pnEl)
      pnEl = pnNode->NextSiblingElement(name);
    return PersistPtr(new XMLPersist(pnEl));
  }

  const char* XMLPersist::XmlRead()const
  {
    return pnNode->GetText();
  }

  ///Look first to see if there is a child element of this name.
  ///Look second for an attribute of this name.
  ///If either found, return its value.
  ///If MustBeThere is true(the default) inserts from defaults.xml and returns its value.
  ///Otherwise returns NULL;
  ///If name is empty, returns NULL if there are no children and an empty string if there are.
  const char* XMLPersist::XmlReadValue(const std::string& name, bool MustBeThere, const char* pluginTypeID)
  {
    const char* ptext=NULL;
    TiXmlElement* pnEl;
    if(name.empty())
      return pnNode->NoChildren() ? NULL : "";

    //Look to see if there is a child element of this name and, if so, return its value
    pnEl = pnNode->FirstChildElement(name);
    if(pnEl)
    {
      ptext = pnEl->GetText();
      if(!ptext)
        return ""; //element exists but is empty
    }
    else
      ptext = pnNode->Attribute(name.c_str());

    if(!ptext && MustBeThere && !name.empty() && InsertDefault(name, pnNode->Value(), pluginTypeID))
      return XmlReadValue(name, MustBeThere); //Try again (recursively). Should succeed.
    return ptext;
  }

  double XMLPersist::XmlReadDouble(const std::string& name, bool MustBeThere)
  {
    double val=NaN; //was 0.0;
    const char* ptxt = XmlReadValue(name, MustBeThere);
    if(ptxt)
    {
      stringstream ss(ptxt);
      ss >> val;
    }
    return val;
  }

  int XMLPersist::XmlReadInteger(const std::string& name, bool MustBeThere)
  {
    int val(0);
    const char* ptxt = XmlReadValue(name, MustBeThere);
    if(ptxt)
    {
      stringstream ss(ptxt);
      ss >> val;
    }
    return val;
  }

  /** 
  Returns the effective content of an CML <property> element
  Looks for child elements pnList of the form:
  <property dictRef="name">
  <scalar> content </scalar>
  </property>
  or alternatively:
  <property title="name">
  <scalar> content </scalar>
  </property>
  In the second case, only  the part of name after any colon (the "localname") is used
  The property can have <array>, <string>, or anything, in place of <scalar>
  The <property> can be the child of a <propertyList> 
  Returns NULL if the appropriate property is not found or if it has no content.
  **/

  const char* XMLPersist::XmlReadProperty(const string& name, bool MustBeThere)
  {
    TiXmlElement* pnPropList = pnNode->FirstChildElement("propertyList");
    if(!pnPropList)
      pnPropList = pnNode; //Be forgiving; we can get by without a propertyList element
    TiXmlElement* pnProp = pnPropList->FirstChildElement("property");

    while(pnProp)
    {
      size_t pos=0;
      const char* pAtt = pnProp->Attribute("dictRef");
      if(!pAtt)
      {
        pAtt = pnProp->Attribute("title");
        //The position of the start of the localname(the bit after the colon)
        pos=name.find(':') + 1; //ok even if there is no colon (if npos=-1)
      }
      if(pAtt && name.compare(pos, name.size()-pos, pAtt)==0)
      {
        TiXmlElement* pnChild = pnProp->FirstChildElement(); //could be <array> or <scalar> or <string>
        if(pnChild)
        {
          const char * ptxt = pnChild->GetText();
          if(ptxt)
            return ptxt;
        }
        cerr << "Malformed property: " << name << endl;
        return NULL;
      }
      pnProp = pnProp->NextSiblingElement();
    }

    // not found
    if(MustBeThere  && InsertDefault("property", name))
      return XmlReadProperty(name, MustBeThere) ; //Try again (recursively). Should succeed.

    return NULL;
  }

  PersistPtr XMLPersist::XmlMoveToProperty(const string& name, bool ToNextProperty)
  {
    // Set ToNextProperty false to search from a parent of <property> or
    // set true to call from <property> or <scalar> of a previous property.
    TiXmlElement* pnProp;
    if(ToNextProperty)
    {
      pnProp = (pnNode->ValueStr()=="property") ? pnNode : pnNode->Parent()->ToElement();
      pnProp = pnProp->NextSiblingElement();
    }
    else
    {
      TiXmlElement* pnPropList = pnNode->FirstChildElement("propertyList");
      if(!pnPropList)
        pnPropList = pnNode; //Be forgiving; we can get by without a propertyList element
      pnProp = pnPropList->FirstChildElement("property");
    }

    while(pnProp)
    {
      size_t pos=0;
      const char* pAtt = pnProp->Attribute("dictRef");
      if(!pAtt)
      {
        pAtt = pnProp->Attribute("title");
        //The position of the start of the localname(the bit after the colon)
        pos=name.find(':') + 1; //ok even if there is no colon (if npos=-1)
      }
      if(pAtt && name.compare(pos, name.size()-pos, pAtt)==0)
        break;
      pnProp = pnProp->NextSiblingElement();
    }
    if(!pnProp)
      return PersistPtr(NULL);

    TiXmlElement* pnChild = pnProp->FirstChildElement(); //returns <array> or <scalar> or <string>
    return PersistPtr(new XMLPersist(pnChild)); 
  }

  double XMLPersist::XmlReadPropertyDouble(const std::string& name, bool MustBeThere)
  {
    double val=NaN; //was 0.0;
    const char* ptxt = XmlReadProperty(name, MustBeThere);
    if(ptxt)
    {
      stringstream ss(ptxt);
      ss >> val;
    }
    return val;
  }

  int XMLPersist::XmlReadPropertyInteger(const std::string& name, bool MustBeThere)
  {
    int val(0);
    const char* ptxt = XmlReadProperty(name, MustBeThere);
    if(ptxt)
    {
      stringstream ss(ptxt);
      ss >> val;
    }
    return val;
  }


  /// Returns the attName attribute of an CML <property> element
  /// See XMLPersist::XmlReadProperty for details
  const char* XMLPersist::XmlReadPropertyAttribute(const string& name, const string& attName, bool MustBeThere)
  {
    TiXmlElement* pnProp = pnNode->FirstChildElement("property");
    while(pnProp)
    {
      size_t pos=0;
      const char* pAtt = pnProp->Attribute("dictRef");
      if(!pAtt)
      {
        pAtt = pnProp->Attribute("title");
        //The position of the start of the localname(the bit after the colon)
        pos=name.find(':') + 1; //ok even if there is no colon (if npos=-1)
      }
      if(pAtt && name.compare(pos, name.size()-pos, pAtt)==0)
      {
        TiXmlElement* pnChild = pnProp->FirstChildElement(); //could be <array> or <scalar> or <string>
        if(pnChild){
          const char* attrtxt = pnChild->Attribute(attName.c_str());
          if (!attrtxt && MustBeThere) //added 17Feb2015
            break;
          //if(attrtxt) Changed CM 13/3/2013 do NOT look on next sibling property if attribute is not present
            return attrtxt;//return if the attribute is on this element; otherwise look on next sibling element
        }
      }
      pnProp = pnProp->NextSiblingElement();
    }

    // element/attribute combination not found
    if(pDefaults && MustBeThere  && InsertDefault("property", name))
    {
      //Try again (recursively). Should succeed.
      const char* ret = XmlReadPropertyAttribute(name, attName); 
      assert(ret);
      return ret;
    }  
    return NULL;
  }

  /// Returns true if datatext associated with name is "1" or "true" or "yes" or nothing;
  //  returns false if datatext is something else or if element is not found.
  bool XMLPersist::XmlReadBoolean( const std::string& name)
  {
    const char* txt = XmlReadValue(name, false);
    if(txt)
    {
      string s(txt);
      return s.empty() || s=="1" || s=="yes" || s=="true";
    }
    return false;
  }

  /// Writes a value to the element
  void XMLPersist::XmlWrite(const std::string& value)
  {
    TiXmlNode* val = pnNode->FirstChild();
    if(!val)
      val = pnNode->InsertEndChild(TiXmlText(value));
    val->SetValue(value);
  }

  PersistPtr XMLPersist::XmlWriteValueElement(const std::string& name,
    const double datum, const int precision, const bool fixedOnly)
  {
    ostringstream sstrdatum ;
    if(precision>=0)
      sstrdatum.precision(precision);
    if(fixedOnly) //no scientific format
      sstrdatum.setf(ios::fixed,ios::floatfield);
    sstrdatum << datum ;


    TiXmlElement* item = new TiXmlElement(name);
    pnNode->LinkEndChild(item);
    item->LinkEndChild(new TiXmlText(sstrdatum.str()));
    return PersistPtr(new XMLPersist(item));
  }

  /// Inserts into XML document a new element  containing a string
  PersistPtr XMLPersist::XmlWriteValueElement(const std::string& name, const std::string& value, bool cdatamode)
  {
    TiXmlElement* item = new TiXmlElement( name );
    pnNode->LinkEndChild(item);
    TiXmlNode* pn = item->LinkEndChild(new TiXmlText(value));
    pn->ToText()->SetCDATA(cdatamode);
    return PersistPtr(new XMLPersist(item));
  }

  PersistPtr XMLPersist::XmlWriteElement(const std::string& name)
  {
    TiXmlElement* item = new TiXmlElement( name );
    pnNode->LinkEndChild(item);
    return PersistPtr(new XMLPersist(item));
  }

  void XMLPersist::XmlWriteAttribute(const std::string& name, const std::string& value)
  {
    pnNode->SetAttribute(name, value);
  }

  void XMLPersist::XmlWriteAttribute(const std::string& name,
    const double datum, const int precision, const bool fixedOnly)
  {
    ostringstream sstrdatum ;
    if(precision>=0)
      sstrdatum.precision(precision);
    if(fixedOnly) //no scientific format
      sstrdatum.setf(ios::fixed,ios::floatfield);
    sstrdatum << datum ;

    pnNode->SetAttribute(name, sstrdatum.str());
  }

  //e.g. <metadata name="dc:source" content="LibraryMols.xml" timestamp="20080705_104810" />
  PersistPtr XMLPersist::XmlWriteMetadata(const std::string& name, const std::string& content)
  {
    TiXmlElement* pn = pnNode->NextSiblingElement("metadataList");
    if(pn)
      pnNode = pn;
    TiXmlElement* item = new TiXmlElement("metadata");
    pnNode->LinkEndChild(item);
    TimeCount events;
    std::string timeString;
    item->SetAttribute("name",name);
    item->SetAttribute("content",content);
    item->SetAttribute("timestamp",events.setTimeStamp(timeString));
    return PersistPtr(new XMLPersist(item));
  }

  PersistPtr XMLPersist::XmlWriteMainElement(
    const std::string& name, const std::string& comment, bool replaceExisting)
  {
    // Delete any existing element of the same name, unless explicitly asked not to.
    if(replaceExisting)
    {
      TiXmlNode* pnExisting = pnNode->FirstChild(name);
      if(pnExisting)
        pnNode->RemoveChild(pnExisting);
    }

    // Add new element with specified name
    TiXmlElement* pnel = new TiXmlElement( name );
    pnNode->LinkEndChild(pnel);

    //No timestamp or comment if comment is empty
    if(comment.size())
    {
      //----------------------------------------
      TimeCount events;
      std::string thisEvent, timeString;
      pnel->SetAttribute("calculated", events.setTimeStamp(timeString));
      //----------------------------------------

      //Add explanatory comment in description element
      TiXmlElement* pncom = new TiXmlElement("me:description");
      pnel->LinkEndChild( pncom );
      pncom->LinkEndChild(new TiXmlText(comment));
    }
    return PersistPtr(new XMLPersist(pnel));
  }

  ///Insert into XML document a new property element
  /**If the parameter units is not empty a timestamp and a units attribute are added. Like:

  <property dictRef="me:ZPE">
  <scalar calculated="20081122_230151" units="kJ/mol">139.5</scalar>
  </property>

  Returns a PersistPtr to the <scalar> element (so that more attributes can be added).
  **/
  PersistPtr XMLPersist::XmlWriteProperty( const std::string& name, 
    const std::string& content, const std::string& units)
  {
    TiXmlElement* pnprop = new TiXmlElement( "property" );
    pnNode->LinkEndChild(pnprop);
    pnprop->SetAttribute("dictRef", name);
    TiXmlElement* pnscal = new TiXmlElement( "scalar" );
    pnprop->LinkEndChild(pnscal);
    pnscal->LinkEndChild(new TiXmlText(content));
    if(!units.empty())
    {
      TimeCount events;
      std::string thisEvent, timeString;
      pnscal->SetAttribute("calculated", events.setTimeStamp(timeString));
      pnscal->SetAttribute("units", units);
    }
    return PersistPtr(new XMLPersist(pnscal));
  }

  PersistPtr XMLPersist::XmlCopy(PersistPtr ppToBeCopied, PersistPtr ppToBeReplaced)
  {
    if(ppToBeReplaced)
    {
      XMLPersist* pxPRep = dynamic_cast<XMLPersist*> (ppToBeReplaced.get());
      if(!pnNode->RemoveChild(pxPRep->pnNode))
        return PersistPtr(NULL); //deletion failure
    }

    TiXmlNode* pnBefore = pnNode->FirstChild();
    if(!pnBefore) return false;
    XMLPersist* pxP = dynamic_cast<XMLPersist*> (ppToBeCopied.get());
    return PersistPtr( new XMLPersist(pnNode->InsertBeforeChild(pnBefore, *pxP->pnNode)->ToElement()));
  }

  bool XMLPersist::XmlSaveFile(const std::string& outfilename)
  {
    if(outfilename.empty())
      return pnNode->GetDocument()->SaveFile(stdout);
    else
      return pnNode->GetDocument()->SaveFile(outfilename);
  }

  bool XMLPersist::InsertDefault(const string& elName, const string& dictRefName, const char* pluginTypeID)
  {
    bool isProperty = elName == "property";
    string name = !isProperty ? elName : dictRefName;

    //Find the default element in defaults.xml
    if (!pDefaults)
    {
      //Open defaults.xml (once only)
      static TiXmlDocument def;
      if (def.LoadFile(m_defaultsfilename))
        pDefaults = &def;
      else
      {
        //Should really terminate, but it's not easy to do from here
        string s("Could not open the defaults file: " + m_defaultsfilename + xtraerr);
        throw (std::runtime_error(s));
      }
    }

    string possibles;
    const char* pattr;
    if (pluginTypeID)
      possibles = TopPlugin::List(pluginTypeID, TopPlugin::comma);

    TiXmlElement* pnDefProp = pDefaults->RootElement()->FirstChildElement(elName);
    // look for a matching element
    if (!isProperty)
    {
      //TiXmlElement* pnDefProp = pDefaults->RootElement()->FirstChildElement(elName);
      if (!isProperty && pnDefProp)
      {
        TiXmlNode* reInserted;
        if (reInserted = pnNode->InsertEndChild(*pnDefProp))
        {
          pattr = pnDefProp->Attribute("default");
          string attrtext = pattr;
          if (attrtext == "true")
            cinfo << "The default value of " << name << " was used." << endl;
          else
          {
            if (pluginTypeID) //write back default attribute with the possible values
              reInserted->ToElement()->SetAttribute("default", attrtext + possibles);
            cerr << "No value of " << name << " was supplied and the default value "
              << attrtext << possibles << endl;
          }
          return true;
        }
      }
      //element not found
      //look for an element dictRefName with attribute elName
      pnDefProp = pDefaults->RootElement()->FirstChildElement(dictRefName);
      while (pnDefProp)
      {
        const char* txt = pnDefProp->Attribute(name.c_str());
        if (txt)
        {
          pnNode->SetAttribute(name, txt);
          pattr = pnDefProp->Attribute("default"); 
          string attrtext(pattr);
          pnNode->SetAttribute("default", attrtext);
          if (attrtext == "true")
            cinfo << "The default value of the attribute " << name << " was used." << endl;
          else
            cerr << "No value of " << name << " was supplied and the default value " << attrtext << endl;
          return true;
        }
        pnDefProp = pnDefProp->NextSiblingElement(dictRefName);
      }
    }
    // look for property
    while (pnDefProp)
    {
      const char* dictRefText = pnDefProp->Attribute("dictRef");
      if (dictRefText && dictRefName == dictRefText)
      {
        pattr = pnDefProp->Attribute("default");
        if (pattr)
        {
          //default on <property> so copy whole element
          pnNode->InsertEndChild(*pnDefProp);
          return true;
        }
        //There is no default attribute on <property> element in defaults.xml,
        //so copy only the attributes on the child element (probably scalar)

        //Find property in datatafile
        TiXmlElement* pnScalar =
          dynamic_cast<XMLPersist*>(XmlMoveToProperty(dictRefName).get())->pnNode;

        TiXmlAttribute* pAttrib = pnDefProp->FirstChildElement()->FirstAttribute();
        while (pAttrib)
        {
          pnScalar->SetAttribute(pAttrib->Name(), pAttrib->Value());
          pAttrib = pAttrib->Next();
        }
        cinfo << "The default attribute(s) were added to the the property " << dictRefName << endl;
        return true;
      }
      pnDefProp = pnDefProp->NextSiblingElement("property");
    }
    cerr << "No element, attribute or dictRef of a property, " << name << ", was found in defaults.xml or it could not be copied."
      << "\nTHIS NEEDS TO BE CORRECTED.\n";
    return false;

  }

}//namespacer mesmer
