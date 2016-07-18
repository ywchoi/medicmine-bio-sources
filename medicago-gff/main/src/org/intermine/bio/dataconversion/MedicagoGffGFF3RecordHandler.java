package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2002-2011 FlyMine
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */

import java.util.Map;
import java.util.List;
import java.util.Arrays;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.apache.commons.lang.StringUtils;

import org.intermine.bio.io.gff3.GFF3Record;
import org.intermine.metadata.Model;
import org.intermine.metadata.StringUtil;
import org.intermine.xml.full.Item;

/**
 * A converter/retriever for the AipGff dataset via GFF files.
 */

public class MedicagoGffGFF3RecordHandler extends GFF3RecordHandler
{

    private final Map<String, Item> protIdMap = new HashMap<String, Item>();

    /**
     * Create a new MedicagoGffGFF3RecordHandler for the given data model.
     * @param model the model for which items will be created
     */
    public MedicagoGffGFF3RecordHandler (Model model) {
        super(model);
        refsAndCollections.put("MRNA", "gene");
        refsAndCollections.put("Exon", "transcripts");
        refsAndCollections.put("FivePrimeUTR", "mRNAs");
        refsAndCollections.put("ThreePrimeUTR", "mRNAs");
        refsAndCollections.put("TRNA", "gene");
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void process(GFF3Record record) {
        // This method is called for every line of GFF3 file(s) being read.  Features and their
        // locations are already created but not stored so you can make changes here.  Attributes
        // are from the last column of the file are available in a map with the attribute name as
        // the key.   For example:
        //
        Item feature = getFeature();
        //     String symbol = record.getAttributes().get("symbol");
        //     feature.setAttrinte("symbol", symbol);
        //
        // Any new Items created can be stored by calling addItem().  For example:
        //
        //     String geneIdentifier = record.getAttributes().get("gene");
        //     gene = createItem("Gene");
        //     gene.setAttribute("primaryIdentifier", geneIdentifier);
        //     addItem(gene);
        //
        // You should make sure that new Items you create are unique, i.e. by storing in a map by
        // some identifier.
        String clsName = feature.getClassName();

        String regexp = "Gene|MRNA";
        Pattern p = Pattern.compile(regexp);
        Matcher m = p.matcher(clsName);
        if (m.find()) {
            List<String> dbxrefs = record.getDbxrefs();
            if (dbxrefs != null) {
                Iterator<String> dbxrefsIter = dbxrefs.iterator();

                while (dbxrefsIter.hasNext()) {
                    String dbxref = dbxrefsIter.next();

                    List<String> refList = new ArrayList<String>(
                            Arrays.asList(StringUtil.split(dbxref, ",")));
                    for (String ref : refList) {
                        ref = ref.trim();
                        int colonIndex = ref.indexOf(":");
                        if (colonIndex == -1) {
                            throw new RuntimeException("external reference not understood: " + ref);
                        }

                        if (ref.startsWith("locus:")) {
                            String locus_tag = ref.substring(colonIndex + 1);
                            feature.setAttribute("secondaryIdentifier", locus_tag);
                        } else if (ref.startsWith("UniProt:")) {
                            String uniprotAcc = ref.substring(colonIndex + 1);

                            Item proteinItem;
                            if (protIdMap.containsKey(uniprotAcc)) {
                                proteinItem = protIdMap.get(uniprotAcc);
                            } else {
                                proteinItem = converter.createItem("Protein");
                                proteinItem.setAttribute("primaryAccession", uniprotAcc);
                                proteinItem.setReference("organism", getOrganism());
                                addItem(proteinItem);

                                protIdMap.put(uniprotAcc, proteinItem);
                            }
                            feature.setReference("protein", proteinItem);

                        } else {
                            throw new RuntimeException("unknown external reference type: " + ref);
                        }
                    }
                }
            }
        }
    }
}
