package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2002-2016 FlyMine
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */

import java.io.File;
import java.io.IOException;
import java.io.Reader;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.apache.tools.ant.BuildException;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.util.FormattedTextParser;
import org.intermine.xml.full.Item;


/**
 *
 * @author Vivek Krishnakumar
 */
public class MtgiGeneMappingConverter extends BioFileConverter
{
    //
    private static final String TAX_ID = "3880";
    private static final String DATASET_TITLE = "MTGI to Gene mapping";
    private static final String DATA_SOURCE_NAME = "DFCI";

    private static final Logger LOG = Logger.getLogger(MtgiGeneMappingConverter.class);

    private Item organism;
    private Map<String, Item> unstoredItems = new HashMap<String, Item>();
    private Map<String, Item> storedItems = new HashMap<String, Item>();

    /**
     * Constructor
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     */
    public MtgiGeneMappingConverter(ItemWriter writer, Model model) throws ObjectStoreException {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
        createOrganismItem();
    }

    /**
     *
     *
     * {@inheritDoc}
     */
    public void process(Reader reader) throws Exception {
        File currentFile = getCurrentFile();
        if (currentFile.getName().contains("map")) {
            processFile(reader, organism);
        }
    }

    /**
     * Process all rows of the Gene Index mapping files
     *
     * @param reader
     * @throws IOException
     * @throws ObjectStoreException
     */
    private void processFile(Reader reader, Item organism)
        throws IOException, ObjectStoreException {
        Iterator<?> tsvIter;
        try {
            tsvIter = FormattedTextParser.parseTabDelimitedReader(reader);
        } catch (Exception e) {
            throw new BuildException("cannot parse file: " + getCurrentFile(), e);
        }

        while (tsvIter.hasNext()) {
            String[] line = (String[]) tsvIter.next();

            String tcIdentifier = line[0];

            // store Tentative Consensus feature
            Item tcFeature = createFeature(tcIdentifier, "TentativeConsensus", organism, false);

            // identify gene feature type based on locus identifier format
            List<String> geneIdentifiers = Arrays.asList(line[1].split(","));
            for (String geneIdentifier : geneIdentifiers) {
                String geneFeatureType = geneIdentifier.contains("te") ? "TransposableElementGene" : "Gene";
                Item geneFeature = createFeature(geneIdentifier, geneFeatureType, organism, false);

                // store reference to TC feature in associatedTCs collection
                geneFeature.addToCollection("associatedTCs", tcFeature.getIdentifier());

                // store reference to gene feature in associatedGenes collection
                tcFeature.addToCollection("associatedGenes", geneFeature.getIdentifier());

                if (!storedItems.containsKey(geneIdentifier)) {
                    store(geneFeature);
                    storedItems.put(geneIdentifier, geneFeature);
                }
            }

            // store tcAliases from line[2] (column 3 of tsv file) if not `null` or "na"
            if (!line[2].equals("na")) {
                List<String> tcAliases = Arrays.asList(line[2].split(","));
                for (String tcAliasId : tcAliases) {
                    Item synonym = storedItems.get(tcAliasId);
                    if (synonym == null) {
                        synonym = createSynonym(tcFeature, tcAliasId, true);
                        storedItems.put(tcAliasId, synonym);
                    }
                    tcFeature.addToCollection("synonyms", synonym);
                }
            }
            if (!storedItems.containsKey(tcIdentifier)) {
                store(tcFeature);
                storedItems.put(tcIdentifier, tcFeature);
            }
        }
    }

    /**
     * Create and store a BioEntity item on the first time called.
     *
     * @param primaryId the primaryIdentifier
     * @param organism the related organism Item
     * @throws ObjectStoreException
     */
    private Item createFeature(String primaryId, String type, Item organism, boolean store)
        throws ObjectStoreException {
        Item feature = null;

        if (storedItems.containsKey(primaryId)) {
            feature = storedItems.get(primaryId);
        } else if (unstoredItems.containsKey(primaryId)) {
            feature = unstoredItems.get(primaryId);
        } else {
            feature = createItem(type);
            feature.setAttribute("primaryIdentifier", primaryId);
            feature.setReference("organism", organism);
            if (store) {
                store(feature);
                storedItems.put(primaryId, feature);
            } else {
                unstoredItems.put(primaryId, feature);
            }
        }

        return feature;
    }

    /**
     * Create and store a organism item on the first time called.
     *
     * @throws ObjectStoreException os
     */
    protected void createOrganismItem() throws ObjectStoreException {
        organism = createItem("Organism");
        organism.setAttribute("taxonId", TAX_ID);
        store(organism);
    }
}
