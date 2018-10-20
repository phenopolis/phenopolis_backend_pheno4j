from views import *
from lookups import *
   
    
@app.route('/login', methods=['POST'])
def login():
    if not request.is_json: return jsonify({"msg": "Missing JSON in request"}), 400
    username = request.json.get('username', None)
    password = request.json.get('password', None)
    if not username: return jsonify({"msg": "Missing username parameter"}), 400
    if not password: return jsonify({"msg": "Missing password parameter"}), 400
    if username != 'test' or password != 'test': return jsonify({"msg": "Bad username or password"}), 401
    # Identity can be any data that is json serializable
    access_token = create_access_token(identity=username)
    return jsonify(access_token=access_token), 200

@app.route('/rv_sharing/<individual_id>/<thresh>/<allele_freq>/<limit>')
@jwt_required
def rv_sharing(individual_id,thresh,allele_freq,limit):
    #thresh=0.05
    #allele_freq=0.001
    print individual_id
    print float(thresh)
    print float(allele_freq)
    neo=get_db('neo4j')
    q= """ MATCH (k:Person)
    WITH count(k) as numberOfPeople
    MATCH (p:Person {{personId:"{personId}"}})<-[:PRESENT_IN]-(gv:GeneticVariant)
    WHERE (gv.allele_freq < {allele_freq} or gv.hasExac = false)
    WITH size(()<-[:PRESENT_IN]-(gv)) as count , gv, p, numberOfPeople
    WHERE count > 1 
    AND ((count / toFloat(numberOfPeople))  <= {thresh})
    MATCH (gv)-[:PRESENT_IN]->(q:Person)
    WHERE p <> q
    WITH p,q,count(gv) as intersection, numberOfPeople
    ORDER BY intersection DESC limit {limit}
    MATCH (x:Person)<-[:PRESENT_IN]-(v:GeneticVariant)
    WHERE (x.personId = p.personId or x.personId = q.personId)
    AND (v.allele_freq < {allele_freq} or v.hasExac = false)
    AND ((size(()<-[:PRESENT_IN]-(v)) / toFloat(numberOfPeople))  <= {thresh})
    WITH p, q, v, intersection
    RETURN p.personId, q.personId, intersection, size(collect(distinct v)) as unionSum, (round((intersection/toFloat(size(collect(distinct v))))*100.0*10)/10) as PercentShared
    ORDER BY PercentShared DESC;
    """.format(thresh=float(thresh), personId=individual_id,allele_freq=float(allele_freq),limit=int(limit))
    result = neo.run(q)
    get_db('neo4j').close()
    return json.dumps([r.__dict__ for r in result], indent=4)


@app.route('/update_patient_data/<individual>',methods=['POST'])
@jwt_required
def update_patient_data(individual):
    if session['user']=='demo': return 'not permitted'
    print(request.form)
    consanguinity=request.form.getlist('consanguinity_edit[]')[0]
    gender=request.form.getlist('gender_edit[]')[0]
    genes=request.form.getlist('genes[]')
    features=request.form.getlist('feature[]')
    print('INDIVIDUAL',individual)
    print('GENDER',gender)
    print('CONSANGUINITY',consanguinity)
    print('GENES',genes)
    print('FEATURES',features)
    print(individual)
    external_id=individual
    individual=get_db(app.config['DB_NAME_PATIENTS']).patients.find_one({'external_id':external_id})
    print('edit patient gender')
    print(get_db(app.config['DB_NAME_PATIENTS']).patients.update_one({'external_id':external_id},{'$set':{'sex':{'female':'F','male':'M','unknown':'U'}[gender]}}))
    print('edit patient genes')
    individual['genes']=[]
    for g in genes:
        gene=get_db(app.config['DB_NAME']).genes.find_one({'gene_name_upper':g})
        print(gene)
        if gene in [g['gene'] for g in individual['genes']]: continue
        if not gene: continue
        individual['genes'].append({'gene':g, 'status':'candidate'})
    print(individual['genes'])
    print(get_db(app.config['DB_NAME_PATIENTS']).patients.update_one({'external_id':external_id},{'$set':{'genes':individual['genes']}}))
    print('edit patient features')
    individual['features']=[]
    for f in features:
        hpo=get_db(app.config['DB_NAME_HPO']).hpo.find_one({'name':re.compile('^'+f+'$',re.IGNORECASE)})
        if not hpo: continue
        if hpo in [h['label'] for h in individual['features']]: continue
        individual['features'].append({'id':hpo['id'][0], 'label':hpo['name'][0], 'observed':'yes'})
    print(get_db(app.config['DB_NAME_PATIENTS']).patients.update_one({'external_id':external_id},{'$set':{'features':individual['features']}}))
    print(get_db(app.config['DB_NAME_PATIENTS']).patients.update_one({'external_id':external_id},{'$set':{'observed_features':[f for f in individual['features'] if f['observed']=='yes']}}))
    print('edit patient consanguinity')
    individual['family_history']=individual.get('family_history',{})
    if (consanguinity)=='unknown':
        individual['family_history']['consanguinity']=None
    elif consanguinity.lower()=='yes':
        individual['family_history']['consanguinity']=True
    elif consanguinity.lower()=='no':
        individual['family_history']['consanguinity']=False
    print(get_db(app.config['DB_NAME_PATIENTS']).patients.update_one({'external_id':external_id},{'$set':{'family_history':individual['family_history']}}))
    # also trigger refresh of that individual for individuals summary page
    patient=Patient(external_id,get_db(app.config['DB_NAME_PATIENTS']))
    print(patient.consanguinity)
    print(patient.observed_features)
    print(patient.genes)
    print(patient.gender)
    individuals_update([external_id])
    return jsonify({'success': True}), 200


@app.route('/venn_json/<individual>')
@jwt_required
def venn_json(individual):
    feature_venn=get_feature_venn(individual)
    return jsonify(result=feature_venn)


@app.route('/homozygous_variants_json/<individual>')
@jwt_required
def homozgous_variants(individual):
    allele_freq=float(request.args.get('allele_freq',0.001))
    kaviar_AF=float(request.args.get('kaviar_AF',0.001))
    s=""" MATCH
    (p)-[:PersonToObservedTerm]-(t:Term),
    (t)--(g:Gene)--(gv:GeneticVariant)-[:HomVariantToPerson]-(p:Person), 
    (gv)--(tv:TranscriptVariant)
    WHERE p.personId='%s' AND gv.kaviar_AF < %f AND gv.allele_freq < %f
    WITH gv, g, t, tv
    OPTIONAL
    MATCH
    (gv)-[:HetVariantToPerson]-(p2:Person)
    OPTIONAL
    MATCH
    (gv)-[:HomVariantToPerson]-(p3:Person)
    RETURN gv,
    collect(distinct g),
    collect(distinct t),
    collect(distinct tv),
    collect(distinct p2),
    collect(distinct p3)
    """ % (individual,kaviar_AF,allele_freq,)
    print(s)
    with neo4j_driver.session() as neo4j_session: 
        result=neo4j_session.run(s)
    return jsonify(result=[merge_dicts(dict(r[0]),
        {'genes':[dict(x) for x in r[1]]},
        {'terms':[dict(x) for x in r[2]]},
        {'transcript_variants':[dict(x) for x in r[3]]},
        {'het_individuals':[dict(x) for x in r[4]]},
        {'hom_individuals':[dict(x) for x in r[5]]}
        ) for r in result])

@app.route('/compound_het_variants_json/<individual>',methods=['GET','POST'])
@jwt_required
def compound_het_variants(individual):
    kaviar_AF=float(request.args.get('kaviar_AF',0.01))
    allele_freq=float(request.args.get('allele_freq',0.01))
    s="""
    MATCH
    (p)-[:PersonToObservedTerm]-(t:Term),
    (g:Gene)--(gv:GeneticVariant)-[:HetVariantToPerson]-(p:Person)
    WHERE p.personId='%s' AND gv.kaviar_AF<%f and gv.allele_freq < %f
    WITH g, collect(distinct gv) AS cgv
    WHERE length(cgv) > 1
    UNWIND cgv as v
    OPTIONAL
    MATCH
    (v)-[:HetVariantToPerson]-(p2:Person)
    OPTIONAL
    MATCH
    (v)-[:HomVariantToPerson]-(p3:Person)
    RETURN v,
    collect(distinct g),
    collect(distinct p2),
    collect(distinct p3)
    """ % (individual,kaviar_AF,allele_freq)
    print(s)
    with neo4j_driver.session() as neo4j_session:
        result=neo4j_session.run(s)
    return jsonify(result=[ merge_dicts(
        dict(r[0]),
        {'terms':[]},
        {'genes':[dict(x) for x in r[1]]},
        {'transcript_variants':[]},
        {'het_individuals':[dict(x) for x in r[2]]},
        {'hom_individuals':[dict(x) for x in r[3]]}
        ) for r in result])


@app.route('/rare_variants_json/<individual>')
@jwt_required
def rare_variants(individual):
    kaviar_AF=float(request.args.get('kaviar_AF',0.01))
    allele_freq=float(request.args.get('allele_freq',0.01))
    s=""" MATCH
    (p)-[:PersonToObservedTerm]-(t:Term),
    (t)--(g:Gene)--(gv:GeneticVariant)-[:HetVariantToPerson]-(p:Person), 
    (gv)--(tv:TranscriptVariant)
    WHERE p.personId='%s' AND gv.kaviar_AF < %f AND gv.allele_freq < %f
    WITH gv, g, t, tv
    OPTIONAL
    MATCH
    (gv)-[:HetVariantToPerson]-(p2:Person)
    OPTIONAL
    MATCH
    (gv)-[:HomVariantToPerson]-(p3:Person)
    RETURN gv,
    collect(distinct g),
    collect(distinct t),
    collect(distinct tv),
    collect(distinct p2),
    collect(distinct p3)
    """ % (individual,kaviar_AF,allele_freq,)
    print(s)
    result=neo4j_session.run(s)
    return jsonify(result=[ merge_dicts(
        dict(r[0]),
        {'genes':[dict(x) for x in r[1]]},
        {'terms':[dict(x) for x in r[2]]},
        {'transcript_variants':[dict(x) for x in r[3]]},
        {'het_individuals':[dict(x) for x in r[4]]},
        {'hom_individuals':[dict(x) for x in r[5]]}
        ) for r in result])
 
@app.route('/individual_update/<individual>')
@jwt_required
def individual_update(individual):
    print 'UPDATE'
    if request.referrer:
        referrer=request.referrer
        u = urlparse(referrer)
        referrer='%s://%s' % (u.scheme,u.hostname,)
        if u.port: referrer='%s:%s' % (referrer,u.port,)
        return redirect(referrer+'/individual/'+individual)
    else:
        return 'done'

@app.route('/homozygous_individuals_json/<variant_id>')
@jwt_required
def get_homozygous_individuals(variant_id):
    s="""
    MATCH
    (v)-[:HomVariantToPerson]-(p:Person)
    WHERE v.variantId='%s'
    RETURN p
    """ % variant_id
    with neo4j_driver.session() as neo4j_session:
        result=neo4j_session.run(s)
    return jsonify(result=[ merge_dicts(
        dict(r[0])) for r in result])


@app.route('/heterozygous_individuals_json/<variant_id>')
@jwt_required
def get_heterozygous_individuals(variant_id):
    s="""
    MATCH
    (v)-[:HetVariantToPerson]-(p:Person)
    WHERE v.variantId='%s'
    RETURN p
    """ % variant_id
    db_session = neo4j_driver.session()
    result=db_session.run(s)
    return jsonify(result=[ merge_dicts(
        dict(r[0])) for r in result])

@app.route('/individual_json/<individual>')
@jwt_required
def individual_json(individual):
    patient=Patient(individual,patient_db=get_db(app.config['DB_NAME_PATIENTS']))
    return patient.json()

@app.route('/hpo_skat_json/<hpo_id>')
@jwt_required
def hpo_skat_json(hpo_id):
    skat_genes=skat(hpo_id)
    return jsonify( result={ 'individuals':skat_genes }, allow_nan=False )

def get_hpo_individuals(hpo_id):
    db=get_db()
    hpo_db=get_db(app.config['DB_NAME_HPO'])
    patients_db=get_db(app.config['DB_NAME_PATIENTS'])
    patients=lookups.get_hpo_patients(hpo_db,patients_db,hpo_id,cached=True)
    print('num patients', len(patients))
    # candidate genes
    candidate_genes = [p.get('genes',[]) for p in patients]
    # solved genes
    solved_genes = [p.get('solved',[]) for p in patients]
    hpo_db=get_db(app.config['DB_NAME_HPO'])
    def f(p):
        print p['external_id']
        if session['user']=='demo': p['external_id']='hidden'
        del p['_id']
        p['features']=[f for f in p.get('features',[]) if f['observed']=='yes']
        if 'solved' in p:
            if 'gene' in p['solved']:
                p['solved']=[p['solved']['gene']]
            else:
                p['solved']=[]
        else: p['solved']=[]
        if 'genes' in p: p['genes']=[x['gene'] for x in p['genes'] if 'gene' in x]
        else: p['genes']=[]
        p['genes']=list(frozenset(p['genes']+p['solved']))
        p2=db.patients.find_one({'external_id':p['external_id']},{'rare_homozygous_variants_count':1,'rare_compound_hets_count':1, 'rare_variants_count':1,'total_variant_count':1})
        if not p2: return p
        p['rare_homozygous_variants_count']=p2.get('rare_homozygous_variants_count','')
        p['rare_compound_hets_count']=p2.get('rare_compound_hets_count','')
        p['rare_variants_count']=p2.get('rare_variants_count','')
        p['total_variant_count']=p2.get('total_variant_count','')
        solved_patient=db.solved_patients.find_one({'external_id':p['external_id']})
        if solved_patient and session['user']!='demo': p['solved_variants']=solved_patient.get('genes',{})
        return p
    patients=[f(p) for p in patients if 'external_id' in p]
    return patients



@app.route('/hpo_individuals_json/<hpo_id>')
@jwt_required
def hpo_individuals_json(hpo_id):
    patients=get_hpo_individuals(hpo_id)
    return jsonify( result={ 'individuals':patients } )


@app.route('/hpo_individuals_csv/<hpo_id>')
@jwt_required
def hpo_individuals_csv(hpo_id):
    patients=get_hpo_individuals(hpo_id)
    return '\n'.join([','.join([p['external_id'],';'.join([str(g) for g in p['genes']])]) for p in patients if 'external_id' in p])


@app.route('/phenogenon_json/<hpo_id>')
@jwt_required
def phenogenon_json(hpo_id):
    cache = bool(request.args.get('cache',True))
    threshold = float(request.args.get('threshold',0.05))
    print 'PHENOGENON_JSON'
    print cache
    #print(intersect(obs_genes.keys(),lit_genes))
    #print(Counter([rv['HUGO'] for rv in db.patients.find_one({'external_id':p['external_id']},{'rare_variants':1})]['rare_variants']))
    ## only return common variants if there are many individuals
    ##rsession.voidEval('common_variants <- common.variants')
    lit_genes=[]
    omim_genes=[]
    recessive_genes=[]
    dominant_genes=[]
    phenogenon(hpo_id, lit_genes, omim_genes, recessive_genes, dominant_genes,cache)
    print(len(lit_genes))
    print(len(omim_genes))
    print(len(recessive_genes))
    print(len(dominant_genes))
    true_positives=len([g for g in lit_genes if 'unrelated_recessive_p_val' in g['phenogenon']['hom_comp'] and g['phenogenon']['hom_comp']['unrelated_recessive_p_val']<=threshold])
    false_negatives=len([g for g in lit_genes if 'unrelated_recessive_p_val' in g['phenogenon']['hom_comp'] and g['phenogenon']['hom_comp']['unrelated_recessive_p_val']>threshold])
    false_positives=len([g for g in recessive_genes if not g['known'] and g['p_val']<=threshold])
    true_negatives=len([g for g in recessive_genes if not g['known'] and g['p_val']>threshold])
    # can have zero denominator sometimes
    #{'TPR':float(true_positives)/float(true_positives+false_negatives),'FPR':float(false_positives)/float(false_positives+true_negatives)}
    return jsonify( result={
        'performance':{'TP':true_positives,'FN':false_negatives,'FP':false_positives,'TN':true_negatives},
                    'lit_genes':lit_genes,
                    'omim_genes':omim_genes,
                    'recessive_genes':recessive_genes,
                    'dominant_genes':dominant_genes,
                    } )


@app.route('/phenogenon_recessive_csv/<hpo_id>')
@jwt_required
def phenogenon_recessive_csv(hpo_id):
    lit_genes=[]
    omim_genes=[]
    recessive_genes=[]
    dominant_genes=[]
    phenogenon(hpo_id, lit_genes, omim_genes, recessive_genes, dominant_genes)
    print(len(lit_genes))
    print(len(omim_genes))
    print(len(recessive_genes))
    print(len(dominant_genes))
    text=','.join([k for k in recessive_genes[0].keys()])+'\n'
    for g in recessive_genes:
        text+=','.join([str(g[k]) for k in recessive_genes[0].keys()])+'\n'
    return text


@app.route('/phenogenon_dominant_csv/<hpo_id>')
@jwt_required
def phenogenon_dominant_csv(hpo_id):
    lit_genes=[]
    omim_genes=[]
    recessive_genes=[]
    dominant_genes=[]
    phenogenon(hpo_id, lit_genes, omim_genes, recessive_genes, dominant_genes)
    print(len(lit_genes))
    print(len(omim_genes))
    print(len(recessive_genes))
    print(len(dominant_genes))
    text=','.join([k for k in dominant_genes[0].keys()])+'\n'
    for g in dominant_genes:
        text+=','.join([str(g[k]) for k in dominant_genes[0].keys()])+'\n'
    return text

@app.route('/phenogenon_literature_csv/<hpo_id>')
@jwt_required
def phenogenon_literature_csv(hpo_id):
    lit_genes=[]
    omim_genes=[]
    recessive_genes=[]
    dominant_genes=[]
    phenogenon(hpo_id, lit_genes, omim_genes, recessive_genes, dominant_genes)
    print(len(lit_genes))
    print(len(omim_genes))
    print(len(recessive_genes))
    print(len(dominant_genes))
    names=['gene_name','phenogenon.dominant_pvalue','phenogenon.recessive_pvalue']
    text=','.join([k for k in names])+'\n'
    for g in lit_genes:
        gene_name=g['gene_name']
        dominant_pvalue=str(g['phenogenon']['het'].get('unrelated_dominant_all_p_val',None))
        recessive_pvalue=str(g['phenogenon']['hom_comp'].get('unrelated_recessive_p_val',None))
        print(gene_name)
        print(dominant_pvalue)
        print(recessive_pvalue)
        text+=','.join([gene_name,dominant_pvalue,recessive_pvalue])+'\n'
    return text


@app.route('/hpo_json/<hpo_id>')
#@auth.login_required
@jwt_required
def hpo_json(hpo_id):
    db=get_db()
    hpo_db=get_db(app.config['DB_NAME_HPO'])
    patients_db=get_db(app.config['DB_NAME_PATIENTS'])
    #patients=[p for p in patients_db.patients.find( { 'features': {'$elemMatch':{'id':str(hpo_id)}} } )]
    print(hpo_id)
    if not hpo_id.startswith('HP:'):
        hpo_term=hpo_db.hpo.find_one({'name':hpo_id})
        hpo_id=hpo_term['id'][0]
    print(hpo_id)
    hpo_term=hpo_db.hpo.find_one({'id':hpo_id})
    hpo_name=hpo_term['name'][0]
    if 'is_a' in hpo_term:
        parents=[ pid for pid in hpo_term['is_a'] ]
    else:
        parents=[]
    print('HPO ANCESTORS')
    hpo_ancestors=lookups.get_hpo_ancestors(hpo_db,hpo_id)
    #print(lookups.get_hpo_ancestors_array(hpo_db,hpo_id))
    print(hpo_ancestors)
    print(len(hpo_ancestors))
    print([h['name'] for h in hpo_ancestors])
    #hpo_ancestors=dict((h['id'][0],h['name'][0]) for h in hpo_ancestors)
    hpo_ancestors=[{'hpo_id':h['id'][0],'hpo_name':h['name'][0]} for h in hpo_ancestors]
    #print(len([v['VARIANT_ID'] for v in db.variants.find({'HET' : { '$in': patient_ids }})]))
    #print(len([v['VARIANT_ID'] for v in db.variants.find({'HOM' : { '$in': patient_ids }})]))
    #r=patients_db.hpo.find_one({'hp_id':hpo_id})
    #if r: external_ids=r['external_ids']
    #else: external_ids=[]
    genes=[lookups.get_gene_by_name(db, r['Gene-Name']) for r in hpo_db.hpo_gene.find({'HPO-ID':hpo_id})]
    print('num genes', len(genes))
    #for r in hpo_db.hpo_pubmed.find({'hpoid':hpo_id}): print(r)
    #pmids=[r['pmid'] for r in hpo_db.hpo_pubmed.find({'hpoid':hpo_id})]
    patients=lookups.get_hpo_patients(hpo_db,patients_db,hpo_id)
    print('num patients', len(patients))
    #return jsonify(result={'hpo_id':hpo_id,'hpo_name':hpo_name,'individuals':[str(p['external_id']) for p in patients],'genes':genes})
    return jsonify(result={'hpo_id':hpo_id,'hpo_name':hpo_name,'individuals':[str(p['external_id']) for p in patients],'hpo_ancestors':hpo_ancestors, 'parents':parents})

@app.route('/gene_json/<gene_id>',methods=['GET','POST'])
@jwt_required
def gene_json(gene_id):
    # if gene not ensembl id then translate to
    db=get_db()
    hpo_db=get_db(app.config['DB_NAME_HPO'])
    patient_db=get_db(app.config['DB_NAME_PATIENTS'])
    hpo=request.args.get('hpo')
    if not gene_id.startswith('ENSG'): gene_id = lookups.get_gene_by_name(get_db(), gene_id)['gene_id']
    gene=db.genes.find_one({'gene_id':gene_id})
    del gene['_id']
    variants=db.variants.find({'genes':gene_id})
    return json.dumps(gene)

@app.route('/gene_phenogenon_json/<gene_id>',methods=['GET','POST'])
@jwt_required
def gene_phenogenon_json(gene_id):
    # if gene not ensembl id then translate to
    db=get_db()
    if not gene_id.startswith('ENSG'): gene_id = lookups.get_gene_by_name(get_db(), gene_id)['gene_id']
    gene=db.gene_hpo_new.find_one({'gene_id':gene_id})
    del gene['_id']
    return json.dumps(gene)

@app.route('/variant_json/<variant_str>')
@jwt_required
def variant_json(variant_str):
    variant=None
    if session['user'] == 'demo':
        variant.__dict__['wt_samples']=[]
        variant.__dict__['het_samples']=[]
        variant.__dict__['hom_samples']=[]
    return jsonify(result=variant.__dict__)

@app.route('/set_variant_causal/<individual>/<variant_str>')
def set_variant_causal(individual, variant_str):
    print individual, variant_str
    db=get_db()
    #get_db().patients.update({'patient_id':individual},{'$addToSet':{'causal_variants':variant_str}})
    var=db.variants.find_one({'variant_id':variant_str})
    gene_id=var['genes'][0]
    gene_name=db.genes.find_one({'gene_id':gene_id})['gene_name_upper']
    print 'GENE_NAME', gene_name
    p=get_db('DB_NAME_PATIENTS').patients.find_one({'external_id':individual})
    get_db('DB_NAME_PATIENTS').patients.update_one({'external_id':individual},{'$set':{'genes': p.get('genes',[])+[{'gene':gene_name}]}})
    print get_db(app.config['DB_NAME_PATIENTS']).patients.update({'external_id':individual},{'$set':p},w=0)
    p=db.patients.find_one({'external_id':individual})
    p['causal_variants']=list(frozenset(p.get('causal_variants',[])+[variant_str]))
    db.patients.update({'external_id':individual},{'$set':{'causal_variants':p['causal_variants']}},w=0)
    if request.referrer:
        referrer=request.referrer
        u = urlparse(referrer)
        referrer='%s://%s' % (u.scheme,u.hostname,)
        if u.port: referrer='%s:%s' % (referrer,u.port,)
    return redirect(referrer+'/individual/'+individual)

@app.route('/unset_variant_causal/<individual>/<variant_str>')
def unset_variant_causal(individual, variant_str):
    print individual, variant_str
    db=get_db()
    p=db.patients.find_one({'external_id':individual})
    if 'causal_variants' in p and not p['causal_variants']: p['causal_variants']=[]
    if variant_str in p.get('causal_variants',[]):
        p['causal_variants']=p['causal_variants'].remove(variant_str)
    db.patients.update({'external_id':individual},{'$set':{'causal_variants':p['causal_variants']}},w=0)
    p2=get_db('DB_NAME_PATIENTS').patients.find_one({'external_id':individual})
    p2['genes']=[]
    for var in p['causal_variants']:
        var=db.variants.find_one({'variant_id':var})
        gene_id=var['genes'][0]
        gene_name=db.genes.find_one({'gene_id':gene_id})['gene_name_upper']
        print 'GENE_NAME', gene_name
        p2['genes']=list(frozenset(p2.get('genes',[])+[{'gene':gene_name}]))
    print get_db(app.config['DB_NAME_PATIENTS']).patients.update({'external_id':individual},{'$set':p2},w=0)
    if request.referrer:
        referrer=request.referrer
        u = urlparse(referrer)
        referrer='%s://%s' % (u.scheme,u.hostname,)
        if u.port: referrer='%s:%s' % (referrer,u.port,)
    return redirect(referrer+'/individual/'+individual)

@app.route('/set_variant_status/<individual>/<variant_str>/<status>')
def set_variant_status(individual, variant_str, status):
    print individual, variant_str, status
    db=get_db()
    #print get_db().patients.update({'patient_id':individual},{'$addToSet':{'variant_status':{variant_str:status}}})
    rare_variants=db.patients.find_one({'external_id':individual},{'rare_variants':1})['rare_variants']
    for rv in rare_variants:
        if rv['variant_id']==variant_str:
            rv['status']=status
    print db.patients.update({'external_id':individual},{'$set':{'rare_variants':rare_variants}})
    return status

@app.route('/private_variants/<individual>')
def private_variants(individual):
    pv=[]
    return jsonify(result=pv)

@app.route('/common_private_variants/<individual>/<individual2>')
def common_private_variants(individual,individual2):
    pv=[]
    return jsonify(result=pv)

@app.route('/common_rare_variants/<individual>/<individual2>/<AC>')
def common_rare_variants(individual,individual2,AC=1):
    pv=[]
    return jsonify(result=pv)

@app.route('/transcript_json/<transcript_id>')
@jwt_required
def transcript_json(transcript_id):
    db = get_db()
    def f(v):
        del v['_id']
        if session['user']=='demo':
            del v['het_samples']
            del v['hom_samples']
            del v['wt_samples']
        return v
    variants=[f(v) for v in db.variants.find({'canonical_transcript':transcript_id})]
    #cache_key = 't-transcript-{}'.format(transcript_id)
    #t = cache.get(cache_key)
    #print 'Rendering %stranscript: %s' % ('' if t is None else 'cached ', transcript_id)
    #if t: return t
    #cache.set(cache_key, t)
    return jsonify(result={'variants':variants,'count':len(variants)})

#@app.route('/individuals/')
def get_individuals(user):
    s="""
    MATCH (u:User {user:'%s'})--(p:Person)-[:PersonToObservedTerm]->(t:Term),
    (p)-[:CandidateGene]-(g:Gene)
    RETURN p.personId as individual,
    p.gender as gender,
    collect(DISTINCT t) as phenotypes,
    p.score as phenotypeScore,
    size((p)<-[:HomVariantToPerson]-()) as hom_count,
    size((p)<-[:HetVariantToPerson]-()) as het_count,
    collect(DISTINCT g.gene_name) as genes;
    """ % user

    with neo4j_driver.session() as db_session: 
        result=db_session.run(s)
    data = []
    for r in result:
        data.append({
            'individual': r['individual'],
            'gender': r['gender'],
            'phenotypes': [dict(x) for x in r['phenotypes']],
            'phenotypeScore': r['phenotypeScore'],
            'hom_count': r['hom_count'],
            'het_count': r['het_count'],
            'genes': [y for y in r['genes']]
        })
    return data


@app.route('/my_patients_json')
@jwt_required
def my_patients_json():
    users_db=get_db(app.config['DB_NAME_USERS'])
    user=users_db.users.find_one({'user':session['user']})
    individuals=get_individuals(user['user'])
    return(jsonify(result=individuals))


# shows each individual, 
# all_individuals
@app.route('/individuals_csv')
def individuals_csv():
    page=int(request.args.get('page',0))
    number=int(request.args.get('number',200))
    hpo_db=get_db(app.config['DB_NAME_HPO'])
    def f(p):
        print p['external_id']
        p['features']=[f for f in p.get('features',[]) if f['observed']=='yes']
        if 'solved' in p:
            if 'gene' in p['solved']:
                p['solved']=[p['solved']['gene']]
            else:
                p['solved']=[]
        else: p['solved']=[]
        if 'genes' in p: p['genes']=[x['gene'] for x in p['genes'] if 'gene' in x]
        else: p['genes']=[]
        p['genes']=list(frozenset(p['genes']+p['solved']))
        p2=get_db().patients.find_one({'external_id':p['external_id']},{'rare_homozygous_variants_count':1,'rare_compound_hets_count':1, 'rare_variants_count':1,'total_variant_count':1})
        if not p2: return p
        p['rare_homozygous_variants_count']=p2.get('rare_homozygous_variants_count','')
        p['rare_compound_hets_count']=p2.get('rare_compound_hets_count','')
        p['rare_variants_count']=p2.get('rare_variants_count','')
        p['total_variant_count']=p2.get('total_variant_count','')
        #p['all_variants_count']=get_db().patients.find_one({'external_id':p['external_id']},{'_id':0,'all_variants_count':1})['all_variants_count']
        #db.cache.find_one({"key" : "%s_blindness,macula,macular,retina,retinal,retinitis,stargardt_" % })
        return p
    conn=PhenotipsClient()
    all_patients=conn.get_patient(session=session).get('patientSummaries',[]) 
    all_eids=[p['eid'] for p in all_patients if p['eid']]
    total=len(all_eids)
    print('TOTAL NUMBER OF PATIENTS',total)
    patients=conn.get_patient(session=session,start=page*number,number=number).get('patientSummaries',[])
    eids=[p['eid'] for p in patients if p['eid']]
    print(eids)
    patients=get_db(app.config['DB_NAME_PATIENTS']).patients.find({'external_id':{'$in':eids}})
    #patients=get_db(app.config['DB_NAME_PATIENTS']).patients.find({'external_id':re.compile('^IRDC')},{'pubmedBatch':0})
    individuals=[f(p) for p in patients if 'external_id' in p]
    # family_history":{"consanguinity":true}
    #if session['user']=='demo': for ind in individuals: ind['external_id']=encrypt(ind['external_id'])
    return '\n'.join([','.join([ind['external_id'],ind['total_variant_count'],ind['rare_variants_count']]) for ind in individuals])


